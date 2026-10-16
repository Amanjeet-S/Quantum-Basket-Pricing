"""Independent reconstruction of the omitted block allowance dictionaries.

Uses independently written residual-envelope formulas. Low-level Arb interval
conversion is the already audited numerical primitive, not the producer's
derivative-envelope or reconstruction function.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,gzip,hashlib,importlib,io,json,sys,time
sys.dont_write_bytecode=True
from flint import arb,ctx


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def reconstruct_one_block(block,norms,bg):
    r=abs(F(block['geometric_loading']))
    M2,M3,M4,D1,D2,D31=[F(norms[name]) for name in ('M2','M3','M4','C2_prime','C2_second','C3_prime')]
    def enclosure(rec):return bg.interval(F(rec['lower_exact_fraction']),F(rec['upper_exact_fraction']))
    def upper(rec):return bg.upper(rec)
    def envelope(kind,ftc=False):
        m={name:enclosure(rec) for name,rec in block['moments'][str(kind)].items()}
        uj=upper((m['U2']*m['J2']).sqrt())
        a=upper((m['U2']*m['U2_J2']).sqrt())
        b=upper((m['U4']*m['J2']).sqrt());u2j=min(a,b)
        if kind in (0,1):
            c={'U_Tj':M2*uj,'r_U2':r*D1*upper(m['U2'])/2,'U2_Tj':M3*u2j/6}
        else:
            c={'Tj2_U_Tjj':M2*(upper(m['J2'])+upper(m['U_JJ'])),
               'r_U_Tj':2*r*D1*uj,'U_Tj2':F(2,3)*M3*upper(m['U_J2']),
               'r2_U2':r*r*D2*upper(m['U2'])/2,'r_U2_Tj':r*D31*u2j/3,
               'U2_Tjj':M3*upper(m['U2_JJ'])/6,'U2_Tj2':M4*upper(m['U2_J2'])/12}
        result={'error_upper':str(sum(c.values(),F(0))),
                'components_upper':{name:str(amount) for name,amount in c.items()},
                'cauchy_variant_a_upper':str(a),'cauchy_variant_b_upper':str(b)}
        if not ftc:result['cauchy_U2_Tj_upper']=str(u2j)
        return result
    return {'block':block['block'],'first_order_FTC':envelope(0,True),
            'centred_first_atoms':envelope(1),'second_Peano':envelope(2)}


def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def streaming_reader_audit(helper):
    sample={'records':[{'index':i,'text':'[ braces } quote " slash \\ newline\n British £ '+str(i),
                       'nested':[None,True,False,-1.25,{'x':[]}]} for i in range(17)],
            'geometry_records':[{'a':[1,2,3]},{'b':'Unicode λ'}],'final':512}
    text=json.dumps(sample,separators=(',',':'),ensure_ascii=True)
    checks=0
    for chunk in (1,2,3,4,7,16,64,1024):
        reader=helper.StreamJSON(io.StringIO(text),chunk=chunk);got={};reader.expect('{');first=True
        while reader.peek()!='}':
            if not first:reader.expect(',')
            name=reader.value();reader.expect(':')
            got[name]=list(reader.entries()) if name in ('records','geometry_records') else reader.value()
            first=False
        reader.expect('}')
        if got!=sample:raise AssertionError(('Streaming reader changed escaped or nested JSON',chunk))
        checks+=1
    # A separate stress test verifies that consumed prefixes are discarded.
    entries=[{'i':i,'payload':'x'*4096} for i in range(80)]
    reader=helper.StreamJSON(io.StringIO(json.dumps(entries,separators=(',',':'))),chunk=1024)
    maximum=0;count=0
    for got in reader.entries():
        maximum=max(maximum,len(reader.buffer));count+=1
        if got!=entries[count-1]:raise AssertionError('Stress reader changed array entry order')
    if count!=len(entries) or maximum>8192:
        raise AssertionError(('Stream buffer grows with whole archive instead of entry size',maximum,count))
    return {'small_chunk_nested_and_escaped_checks':checks,'stress_entries':count,
            'maximum_stress_buffer_characters':maximum,'largest_payload_characters':4096,'status':'passed'}


def audit_omitted_blocks(fundamental,bg):
    manifest=fundamental['reconstruction_manifest'];models={m['model_id']:m for m in MODEL_PROTOCOL['models']}
    raw_geometry=fundamental['geometry_records']
    geom={(models[r['model_id']]['assets'],r['B']):r for r in raw_geometry}
    if len(geom)!=len(raw_geometry):raise AssertionError('Raw geometry has duplicate coverage')
    hashes=manifest['record_semantic_hashes'];records=fundamental['records']
    if len(hashes)!=180 or len(records)!=180:raise AssertionError('All180 payoff records not retained')
    coordinates=blocks=matches=components=0
    for record,expected in zip(records,hashes):
        if (record['model_id'],record['cap'],record['B'])!=(expected['model_id'],expected['cap'],expected['B']):
            raise AssertionError('Record identity or list order changed')
        geometry=geom[(models[record['model_id']]['assets'],record['B'])]
        if len(record['coordinates'])!=len(expected['coordinates']) or len(record['coordinates'])!=len(geometry['coordinates']):
            raise AssertionError('Coordinate counts differ under reconstruction')
        rebuilt={k:v for k,v in record.items() if k!='coordinates'};rebuilt['coordinates']=[]
        for retained,raw,coordinate_hash in zip(record['coordinates'],geometry['coordinates'],expected['coordinates']):
            if 'blocks' in retained:raise AssertionError('Compact record retained an unexpectedly duplicated block list')
            if retained['coordinate']!=raw['coordinate'] or retained['coordinate']!=coordinate_hash['coordinate']:
                raise AssertionError('Coordinate identity or ordering changed')
            omitted=[reconstruct_one_block(b,record['norms'],bg) for b in raw['blocks']]
            if len(omitted)!=coordinate_hash['block_count'] or len(omitted)!=record['B']:
                raise AssertionError('Omitted block count incorrect')
            if canonical_sha(omitted)!=coordinate_hash['canonical_semantic_sha256']:
                raise AssertionError('Independent envelope formulas failed exact block-list reconstruction')
            sums={'0':F(0),'1':F(0),'2':F(0)}
            for block in omitted:
                for kind,field in (('0','first_order_FTC'),('1','centred_first_atoms'),('2','second_Peano')):
                    sums[kind]+=F(block[field]['error_upper']);components+=len(block[field]['components_upper'])
            if {k:str(v) for k,v in sums.items()}!=retained['component_sums_upper']:
                raise AssertionError('Reconstructed block allowances disagree with retained coordinate sums')
            rebuilt['coordinates'].append({**retained,'blocks':omitted})
            coordinates+=1;blocks+=len(omitted);matches+=1
        if canonical_sha(rebuilt)!=expected['canonical_semantic_sha256']:
            raise AssertionError('Complete payoff record did not reconstruct exactly')
        matches+=1
    return {'status':'passed','payoff_records':len(records),'coordinates':coordinates,
            'independently_reconstructed_block_records':blocks,'reconstructed_component_values':components,
            'exact_canonical_semantic_hash_matches':matches,
            'compressed_original_byte_recovery_claim':False,
            'scope':'Own FTC, centred-first and Peano envelope formulas recover every omitted string exactly at original512-bit Arb precision; all aggregate coordinate sums recover exactly.'}


def audit_original_preservation(original,fundamental,helper):
    manifest=fundamental['reconstruction_manifest']
    if sha_file(original)!=manifest['original_archive_sha256'] or original.stat().st_size!=manifest['original_archive_bytes']:
        raise AssertionError('Original compressed archive differs from retained byte checksum')
    seen=[];record_count=geometry_count=field_count=0
    with gzip.open(original,'rt') as inp:
        reader=helper.StreamJSON(inp);reader.expect('{');first=True
        while reader.peek()!='}':
            if not first:reader.expect(',')
            name=reader.value();reader.expect(':');seen.append(name)
            if name=='records':
                for obj in reader.entries():
                    retained=fundamental['records'][record_count];expected=manifest['record_semantic_hashes'][record_count]
                    if canonical_sha(obj)!=expected['canonical_semantic_sha256']:
                        raise AssertionError('Manifest full-record hash does not match actual original')
                    for coordinate,coordinate_hash in zip(obj['coordinates'],expected['coordinates']):
                        omitted=coordinate.pop('blocks')
                        if canonical_sha(omitted)!=coordinate_hash['canonical_semantic_sha256']:
                            raise AssertionError('Manifest omitted-list hash does not match actual original')
                    if obj!=retained:raise AssertionError('A non-omitted payoff record field changed or disappeared')
                    record_count+=1
            elif name=='geometry_records':
                for obj in reader.entries():
                    retained=fundamental['geometry_records'][geometry_count];expected=manifest['geometry_semantic_hashes'][geometry_count]
                    if obj!=retained or canonical_sha(obj)!=expected['canonical_semantic_sha256']:
                        raise AssertionError('Raw geometry evidence changed under compaction')
                    geometry_count+=1
            else:
                obj=reader.value()
                if obj!=fundamental[name] or canonical_sha(obj)!=manifest['preserved_top_field_semantic_hashes'][name]:
                    raise AssertionError(('Original input, mass, timing, counter or other top field changed',name))
                field_count+=1
            first=False
        reader.expect('}')
    if seen!=manifest['original_top_level_key_order']:
        raise AssertionError('Original top-level order or field coverage changed')
    if set(fundamental)!=set(seen)|{'reconstruction_manifest'}:
        raise AssertionError('Compaction contains unexpected or missing top-level fields')
    if record_count!=180 or geometry_count!=len(fundamental['geometry_records']):
        raise AssertionError('Original array coverage incomplete')
    return {'status':'passed','original_records_checked':record_count,'raw_geometry_records_checked':geometry_count,
            'other_top_level_fields_checked':field_count,'original_archive_sha256':manifest['original_archive_sha256'],
            'retention':'Only records[*].coordinates[*].blocks removed; every other original field compared literally.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--compact',type=Path,required=True);ap.add_argument('--original',type=Path);ap.add_argument('--out',type=Path)
    args=ap.parse_args();start_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    sys.path.insert(0,str(args.root));bg=importlib.import_module('boundary_grid_certificate')
    helper=importlib.import_module('block_certificate_archive');helper_sha=sha_file(args.root/'block_certificate_archive.py')
    archives=importlib.import_module('artifact_archive')
    global MODEL_PROTOCOL
    MODEL_PROTOCOL=json.loads((args.root/'results/generalisation_protocol.json').read_text())
    ctx.prec=512
    with archives.open_archive(args.compact) as f:data=json.load(f)
    if data['reconstruction_manifest']['helper_source_sha256']!=helper_sha:
        raise AssertionError('Compact archive was generated with a different helper source')
    parser=streaming_reader_audit(helper)
    result=audit_omitted_blocks(data,bg)
    result['streaming_reader_audit']=parser
    result['original_preservation_audit']=audit_original_preservation(args.original,data,helper) if args.original else {
        'status':'not requested','scope':'Independent reconstruction checked against preserved semantic hashes; raw original byte/field comparison requires --original.'}
    for name,digest in data['source_sha256'].items():
        if sha_file(args.root/name)!=digest:raise AssertionError('Original scientific source has changed')
    if sha_file(args.root/'block_certificate_archive.py')!=helper_sha:
        raise AssertionError('Archive helper changed during independent review')
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=start_sha:
        raise AssertionError('Independent archive auditor changed during execution')
    result['source_sha256']=start_sha;result['helper_source_sha256']=helper_sha;result['compact_sha256']=archives.archive_sha256(args.compact)
    out=args.out or args.root/'results/block_fundamental_independent_audit.json'
    out.write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if not isinstance(v,(list,dict))})


if __name__=='__main__':main()

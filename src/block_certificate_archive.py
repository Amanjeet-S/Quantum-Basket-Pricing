"""Preserve raw block certificate evidence without duplicated payoff strings.

Compaction retains all original fields except coordinates[*].blocks in the
180 payoff records. These lists are reconstructed exactly from retained raw
geometry intervals and six retained rational boundary constants. Canonical
JSON semantic hashes verify reconstruction; compressed original bytes are
not claimed recoverable. The full original archive remains in staging.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json,re,time
from flint import ctx
from artifact_archive import open_archive,archive_sha256
import tensor_transfer as tt
import block_certificate_comparison as bc

ROOT=Path(__file__).resolve().parent
CANON=json.JSONEncoder(sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False)


def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def semantic_sha(obj):
    h=hashlib.sha256()
    for s in CANON.iterencode(obj):h.update(s.encode('utf-8'))
    return h.hexdigest()


class StreamJSON:
    """Sequential top-level reader; at most one nested array entry in memory."""
    tokens=re.compile(r'["\\{}\[\]]')
    def __init__(self,handle,chunk=1<<20):
        self.f=handle;self.chunk=chunk;self.buffer='';self.pos=0;self.eof=False
    def refill(self):
        if self.eof:return False
        if self.pos:
            self.buffer=self.buffer[self.pos:];self.pos=0
        s=self.f.read(self.chunk)
        if not s:self.eof=True;return False
        self.buffer+=s;return True
    def ws(self):
        while True:
            while self.pos<len(self.buffer) and self.buffer[self.pos].isspace():self.pos+=1
            if self.pos<len(self.buffer):return
            if not self.refill():raise EOFError('Unexpected end of JSON')
    def expect(self,ch):
        self.ws()
        if self.buffer[self.pos]!=ch:raise ValueError('Expected '+ch)
        self.pos+=1
    def peek(self):self.ws();return self.buffer[self.pos]
    def value(self):
        self.ws()
        # Discard every previously consumed prefix before a new nested value.
        # Individual records can be large, but preceding records never remain.
        if self.pos:self.buffer=self.buffer[self.pos:];self.pos=0
        start=self.pos;first=self.buffer[start]
        if first not in '{["':
            while True:
                hit=re.search(r'[,}\]\s]',self.buffer[self.pos:])
                if hit:
                    end=self.pos+hit.start();out=json.loads(self.buffer[start:end]);self.pos=end;return out
                # Scalars occur only in small metadata; preserve their start.
                s=self.f.read(self.chunk)
                if not s:raise EOFError('Truncated scalar')
                self.buffer+=s
        scan=start;depth=0;string=False;skip=0
        while True:
            found=self.tokens.search(self.buffer,scan)
            if found is None:
                scan=len(self.buffer)
                s=self.f.read(self.chunk)
                if not s:raise EOFError('Truncated compound JSON value')
                self.buffer+=s;continue
            at=found.start();ch=found.group();scan=at+1
            if at<skip:continue
            if string:
                if ch=='\\':skip=at+2
                elif ch=='"':
                    string=False
                    if first=='"' and depth==0:
                        end=at+1;out=json.loads(self.buffer[start:end]);self.pos=end;return out
            elif ch=='"':string=True
            elif ch in '{[':depth+=1
            elif ch in '}]':
                depth-=1
                if depth==0:
                    end=at+1;out=json.loads(self.buffer[start:end]);self.pos=end;return out
    def entries(self):
        self.expect('[');first=True
        while self.peek()!=']':
            if not first:self.expect(',')
            yield self.value();first=False
        self.expect(']')


def compact(original,fundamental,manifest_path=None):
    ctx.prec=512;start=time.perf_counter();manifest={
        'schema':'qbp-block-evidence-reconstruction-v1',
        'canonical_JSON':'UTF-8 JSON; sort_keys=True, separators=(comma,colon), ensure_ascii=True, allow_nan=False; list order retained; exact strings and parsed numeric values retained',
        'omitted_fields':'Only records[*].coordinates[*].blocks; all raw geometry records, cell/block masses, certificates, aggregate/coordinate component sums, adverse bounds/statuses, norms, timings and counters retained',
        'original_archive_file':Path(original).name,'original_archive_sha256':sha_file(original),
        'original_archive_bytes':Path(original).stat().st_size,'original_top_level_key_order':[],
        'record_semantic_hashes':[],'geometry_semantic_hashes':[],'preserved_top_field_semantic_hashes':{},
        'helper_source_sha256':sha_file(__file__),
        'archive_access_helper_sha256':sha_file(ROOT/'artifact_archive.py'),
        'reconstruction_formula_source_sha256':sha_file(ROOT/'block_certificate_comparison.py'),
        'version_requirement':'Same checked Python/Arb interval conversion and precision512; scientific identity additionally has an independent reconstruction audit',
        'byte_recovery_claim':False,
        'input_certificate_sha256':{str(p.relative_to(ROOT)):sha_file(p) for p in (
          ROOT/'results/boundary_derivative_norms.json',ROOT/'results/gaussian_derivative_norms.json',ROOT/'results/generalisation_protocol.json')}}
    fundamental=Path(fundamental);fundamental.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(original,'rt') as inp,gzip.open(fundamental,'wt') as out:
        reader=StreamJSON(inp);reader.expect('{');out.write('{');first=True
        while reader.peek()!='}':
            if not first:reader.expect(',');out.write(',')
            name=reader.value();reader.expect(':');manifest['original_top_level_key_order'].append(name)
            out.write(json.dumps(name)+':')
            if name in ('records','geometry_records'):
                out.write('[');entryfirst=True
                for obj in reader.entries():
                    if not entryfirst:out.write(',')
                    if name=='records':
                        original_sha=semantic_sha(obj);blocks=[]
                        for coordinate in obj['coordinates']:
                            omitted=coordinate.pop('blocks')
                            blocks.append({'coordinate':coordinate['coordinate'],'block_count':len(omitted),'canonical_semantic_sha256':semantic_sha(omitted)})
                        manifest['record_semantic_hashes'].append({'model_id':obj['model_id'],'cap':obj['cap'],'B':obj['B'],
                          'canonical_semantic_sha256':original_sha,'coordinates':blocks})
                    else:manifest['geometry_semantic_hashes'].append({'model_id':obj['model_id'],'B':obj['B'],'canonical_semantic_sha256':semantic_sha(obj)})
                    json.dump(obj,out,separators=(',',':'));entryfirst=False
                out.write(']')
            else:
                obj=reader.value();manifest['preserved_top_field_semantic_hashes'][name]=semantic_sha(obj)
                json.dump(obj,out,separators=(',',':'))
            first=False
        reader.expect('}');manifest['compaction_seconds_diagnostic']=time.perf_counter()-start
        out.write(',"reconstruction_manifest":');json.dump(manifest,out,separators=(',',':'));out.write('}')
    if manifest_path:
        Path(manifest_path).write_text(json.dumps({**manifest,'fundamental_archive_sha256':sha_file(fundamental),
           'fundamental_archive_file':fundamental.name,'fundamental_archive_bytes':fundamental.stat().st_size},indent=2)+'\n')
    return manifest


def norm_dictionary(record):
    """Recover precisely the six sufficient envelope constants, algebraically."""
    n=record['norms'];M2=F(n['M2']);D1=F(n['C2_prime']);D2=F(n['C2_second'])
    d={'payoff':record['payoff'],'M2_upper':n['M2'],'C_r_upper':{'3':n['M3'],'4':n['M4']}}
    if record['payoff']=='capped':
        d['scalar_boundary_h_first_upper']=str(D1/(2*M2));d['scalar_boundary_h_second_upper']=str(D2/(2*M2))
    elif record['payoff']=='vanilla':
        d['scalar_boundary_density_first_upper']=str(D1/M2);d['scalar_boundary_density_second_upper']=str(D2/M2)
    else:raise ValueError('Unknown retained payoff')
    if min(F(n['M4'])+3*F(n['M3']),D2+2*D1)!=F(n['C3_prime']):raise ValueError('Stored boundary recurrence constant inconsistent')
    return d


def reconstructed_record(record,geometry):
    """Literal frozen envelope functions; no finance-price calculation."""
    ctx.prec=512;norms=norm_dictionary(record);out={k:v for k,v in record.items() if k!='coordinates'};coords=[]
    if len(record['coordinates'])!=len(geometry['coordinates']):raise ValueError('Coordinate count mismatch')
    for retained,raw in zip(record['coordinates'],geometry['coordinates']):
        if retained['coordinate']!=raw['coordinate']:raise ValueError('Coordinate order mismatch')
        blocks=[]
        for block in raw['blocks']:
            envelopes,n=tt.derivative_envelope(norms,block)
            if n!=record['norms']:raise ValueError('Sufficient norm constants changed under reconstruction')
            blocks.append({'block':block['block'],'first_order_FTC':bc.first_only(norms,block,0),
              'centred_first_atoms':envelopes['1'],'second_Peano':envelopes['2']})
        coords.append({**retained,'blocks':blocks})
    out['coordinates']=coords;return out


def verify_and_expand(fundamental,expanded=None,report=None):
    ctx.prec=512;start=time.perf_counter()
    with open_archive(fundamental) as f:data=json.load(f)
    manifest=data['reconstruction_manifest'];stats={'status':'passed','records':0,'coordinates':0,'blocks':0,'semantic_hash_matches':0,
      'original_archive_sha256':manifest['original_archive_sha256'],'compressed_byte_recovery_claim':False}
    if sha_file(ROOT/'block_certificate_comparison.py')!=manifest['reconstruction_formula_source_sha256']:raise ValueError('Frozen producer source differs')
    for field,h in manifest['preserved_top_field_semantic_hashes'].items():
        if semantic_sha(data[field])!=h:raise ValueError('Preserved original field changed: '+field)
    geometry={(r['model_id'].split('_')[0],r['B']):r for r in data['geometry_records']}
    for item in manifest['geometry_semantic_hashes']:
        key=(item['model_id'].split('_')[0],item['B'])
        if semantic_sha(geometry[key])!=item['canonical_semantic_sha256']:raise ValueError('Raw geometry evidence changed')
    def records():
        for r,h in zip(data['records'],manifest['record_semantic_hashes']):
            if (r['model_id'],r['cap'],r['B'])!=(h['model_id'],h['cap'],h['B']):raise ValueError('Record identity/order changed')
            full=reconstructed_record(r,geometry[(r['model_id'].split('_')[0],r['B'])])
            for coordinate,coordinate_hash in zip(full['coordinates'],h['coordinates']):
                blocks=coordinate['blocks']
                if len(blocks)!=coordinate_hash['block_count'] or semantic_sha(blocks)!=coordinate_hash['canonical_semantic_sha256']:
                    raise ValueError('Omitted block fields failed exact semantic reconstruction')
                stats['coordinates']+=1;stats['blocks']+=len(blocks);stats['semantic_hash_matches']+=1
            if semantic_sha(full)!=h['canonical_semantic_sha256']:raise ValueError('Full record failed semantic reconstruction')
            stats['records']+=1;stats['semantic_hash_matches']+=1
            yield full
    if len(data['records'])!=180 or len(manifest['record_semantic_hashes'])!=180:raise ValueError('Original180 outcomes not preserved')
    if expanded:
        with gzip.open(expanded,'wt') as out:
            out.write('{')
            for index,key in enumerate(manifest['original_top_level_key_order']):
                if index:out.write(',')
                out.write(json.dumps(key)+':')
                if key=='records':
                    out.write('[')
                    for k,r in enumerate(records()):
                        if k:out.write(',')
                        json.dump(r,out,separators=(',',':'))
                    out.write(']')
                else:json.dump(data[key],out,separators=(',',':'))
            out.write('}')
    else:
        for _ in records():pass
    stats['verification_seconds_diagnostic']=time.perf_counter()-start
    stats['fundamental_archive_sha256']=archive_sha256(fundamental)
    if report:Path(report).write_text(json.dumps(stats,indent=2)+'\n')
    return stats


def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='mode',required=True)
    cp=sub.add_parser('compact');cp.add_argument('--original',type=Path,required=True);cp.add_argument('--out',type=Path,required=True);cp.add_argument('--manifest',type=Path)
    vp=sub.add_parser('verify');vp.add_argument('--fundamental',type=Path,required=True);vp.add_argument('--expand',type=Path);vp.add_argument('--report',type=Path)
    a=ap.parse_args()
    if a.mode=='compact':r=compact(a.original,a.out,a.manifest);print('Compacted fundamental evidence. Original compressed bytes remain separate.',flush=True)
    else:r=verify_and_expand(a.fundamental,a.expand,a.report);print(json.dumps(r,indent=2),flush=True)

if __name__=='__main__':main()

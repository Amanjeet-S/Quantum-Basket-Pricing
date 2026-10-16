"""Literal retained-gate replay; needs NumPy, without synthesis packages."""
import math


def replay(tape,state):
    import numpy as np
    state=np.asarray(state,dtype=complex).copy();N=len(state);indices=np.arange(N)
    for g in tape:
        kind=g[0]
        if kind=='RY':
            wire=g[1];theta=float(g[2]);a=indices[(indices&(1<<wire))==0];b=a|(1<<wire)
            x,y=state[a].copy(),state[b].copy();c,s=math.cos(theta/2),math.sin(theta/2)
            state[a],state[b]=c*x-s*y,s*x+c*y
        elif kind in ('X','CX','CCX'):
            *controls,target=g[1:];a=indices[(indices&(1<<target))==0]
            for wire in controls:a=a[(a&(1<<wire))!=0]
            b=a|(1<<target);state[a],state[b]=state[b].copy(),state[a].copy()
        elif kind in ('Z','CZ'):
            a=indices
            for wire in g[1:]:a=a[(a&(1<<wire))!=0]
            state[a]*=-1
        elif kind=='H':
            wire=g[1];a=indices[(indices&(1<<wire))==0];b=a|(1<<wire)
            x,y=state[a].copy(),state[b].copy();state[a],state[b]=(x+y)/math.sqrt(2),(x-y)/math.sqrt(2)
        elif kind in ('S','Sdg','T','Tdg'):
            phase={'S':1j,'Sdg':-1j,'T':complex(math.cos(math.pi/4),math.sin(math.pi/4)),'Tdg':complex(math.cos(-math.pi/4),math.sin(-math.pi/4))}[kind]
            state[(indices&(1<<g[1]))!=0]*=phase
        else:raise ValueError(kind)
    return state


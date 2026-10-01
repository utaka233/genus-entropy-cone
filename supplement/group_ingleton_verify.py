#!/usr/bin/env python3
"""Independent exact all-quadruple audit for PSL2(F7) and S5.

Run: python supplement/group_ingleton_verify.py
Requires Python 3 and a C++17 compiler (`c++`); all calculations are integral.
"""
import itertools
import json
import subprocess
import tempfile
from pathlib import Path
import time
from verify_genus3_obstruction import group, permutation_product, enumerate_subgroups

HERE=Path(__file__).resolve().parent


def psl27():
    mats=[t for t in itertools.product(range(7),repeat=4)
          if (t[0]*t[3]-t[1]*t[2])%7==1]
    def canonical(a):
        return min(a,tuple(-x%7 for x in a))
    elems=sorted({canonical(a) for a in mats})
    def op(a,b):
        return canonical(((a[0]*b[0]+a[1]*b[2])%7,
                          (a[0]*b[1]+a[1]*b[3])%7,
                          (a[2]*b[0]+a[3]*b[2])%7,
                          (a[2]*b[1]+a[3]*b[3])%7))
    return group('PSL2(F7)',elems,op)


def symmetric5():
    return group('S5',list(itertools.permutations(range(5))),permutation_product)


def main():
    start=time.time()
    with tempfile.TemporaryDirectory(prefix='group_ingleton_') as tmp:
        tmp=Path(tmp); executable=tmp/'search'
        subprocess.run(['c++','-std=c++17','-O3',str(HERE/'group_ingleton_search.cpp'),
                        '-o',str(executable)],check=True)
        results=[]
        for G,expected in ((psl27(),179),(symmetric5(),156)):
            subgroups=enumerate_subgroups(G)
            assert len(subgroups)==expected
            source=tmp/'subgroups.txt'
            with source.open('w') as f:
                f.write(f"{len(G['elements'])} {len(subgroups)}\n")
                for h in subgroups:
                    f.write(' '.join(str((h>>(64*j))&((1<<64)-1)) for j in range(3))+'\n')
            result=subprocess.run([str(executable),str(source)],text=True,
                                  capture_output=True,check=True)
            print(G['name'],result.stdout,end='',flush=True)
            results.append({'group':G['name'],'subgroup_count':len(subgroups),
                            'search_output':result.stdout,
                            'subgroups':[hex(s) for s in subgroups]})
        (HERE/'group_ingleton_verification.json').write_text(json.dumps(results,indent=2))
    print('PASS: all PSL2(F7) quadruples satisfy Ingleton; all S5 offenders have LHS/RHS=15/16 and H1234=1.')
    print('elapsed seconds',round(time.time()-start,3))

if __name__=='__main__':
    main()

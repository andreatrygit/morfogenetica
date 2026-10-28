import json,sys
d=json.load(open(sys.argv[1]))
W,H=78,36
xs=[p[0] for p in d]; ys=[p[1] for p in d]
lo,hi=-75,75
g=[[' ']*W for _ in range(H)]
for x,y,m in d:
    # x = lungo l'asse del polo (+ = verso il polo), y = trasversale
    cx=int((x-lo)/(hi-lo)*(W-1)); cy=int((y-lo)/(hi-lo)*(H-1))
    if 0<=cx<W and 0<=cy<H: g[H-1-cy][cx]='#' if m else 'o'
print('  asse del polo → destra;  # = cellule costrette')
print('+'+'-'*W+'+')
for r in g: print('|'+''.join(r)+'|')
print('+'+'-'*W+'+')

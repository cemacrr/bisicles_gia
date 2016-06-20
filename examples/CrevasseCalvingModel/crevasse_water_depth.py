
def crevasse_water_depth(x,y,t,thck,topg,*etc):
    
    depth = 0.0
    if (t >= 50.0):
        depth = 25.0
    if (t >= 100.0):
        depth = 20.0
    if (t >= 550.0):
        depth = 0.0

    return depth 

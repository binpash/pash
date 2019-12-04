## Copied from http://code.activestate.com/recipes/577225-union-find/

## This is an interesting way to do a union find because it saves the
## information about the parents in the nodes themselves. I am not
## sure if this can lead to any problem down the line.
##
## Question: Is there any problem that could appear because of this
##           Union-Find implementation?

## TODO: We probably want to do the matching on the ident of the
##       parent and self and not on the object

def MakeSet(x):
     x.parent = x
     x.rank   = 0

def Union(x, y):
     xRoot = Find(x)
     yRoot = Find(y)
     if xRoot.rank > yRoot.rank:
         yRoot.parent = xRoot
     elif xRoot.rank < yRoot.rank:
         xRoot.parent = yRoot
     elif xRoot != yRoot: # Unless x and y are already in same set, merge them
         yRoot.parent = xRoot
         xRoot.rank = xRoot.rank + 1

def Find(x):
     if x.parent == x:
        return x
     else:
        x.parent = Find(x.parent)
        return x.parent

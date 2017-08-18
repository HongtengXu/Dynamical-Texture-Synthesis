function path=PathSelection(W)


[Dis,sp] = dijkstra(sparse(W), 1);
path=[];
loc=sp(end);
while loc~=-1
    path=[loc,path];
    loc=sp(loc);
end
    
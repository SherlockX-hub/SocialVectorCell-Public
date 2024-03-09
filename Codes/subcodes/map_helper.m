function [info,smap,map]=map_helper(x1,y1,frB,dmap)
N=size(frB,2);
% dmap=2;
xbins=0:dmap:70; xbins(1) = [];
[X,Y]=meshgrid(round(xbins/dmap));
% x1r=round(x1/dmap);
% y1r=round(y1/dmap);

x1r=ceil(x1/dmap); %Xiang uses ceil
y1r=ceil(y1/dmap);
gw=1;

occ=cell2mat(arrayfun(@(z) sum((x1r==X(z) & y1r==Y(z))),1:numel(X),'UniformOutput',0)');
socc=general.smooth(reshape(occ,numel(xbins),[]),gw);
socc(occ==0)=nan;

map=cell2mat(arrayfun(@(z) sum(frB(x1r==X(z) & y1r==Y(z),:),1),1:numel(X),'UniformOutput',0)');
smap=arrayfun(@(z) general.smooth(reshape(map(:,z),numel(xbins),[]),gw)./socc,1:N,'UniformOutput',0);
map=map./occ;
% info=arrayfun(@(z) spatial_info_rate(map(:,z)./occ,occ),1:N);
info=arrayfun(@(z) spatial_info_rate(smap{z},socc),1:N);

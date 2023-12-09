% Arc ASC to xyz
% Only works for a raster ASC file converted from ArcMap

clear all
close all

fn = 'grid.asc';
fid = fopen(fn);

% open output file
fidw = fopen(sprintf('%s.xyz',fn(1:end-4)),'wt');
%% read header
S = textscan(fid,'%*s %f',6);
ncols = S{1}(1);
nrows = S{1}(2);
xcor = S{1}(3);
ycor = S{1}(4);
dxy = S{1}(5);
NoData = S{1}(6);

%% read raster data
ffmt = repmat('%f',1,ncols);

for n = nrows:-1:1
   C = textscan(fid,ffmt,1); 
   temp = cell2mat(C);
   
   id = find(temp~=-9999);
   if isempty(id)
       continue
   end
   xcoord = id.*dxy - dxy/2 + xcor;
   ycoord = zeros(1,length(xcoord)) + dxy*n - dxy/2 + ycor;
   zcoord = temp(id);
   
   fprintf(fidw,'%.1f %.1f %.2f\n',[xcoord' ycoord' zcoord']');
   clear C
   
   if mod(n,10) ==0
       if ~isempty(id)
           display(sprintf('%3.2f %% Done - %d points extracted',(nrows-n)/nrows*100,length(id)))
       else
           display(sprintf('%3.2f %% Done',(nrows-n)/nrows*100))
       end
   end

end

fclose('all');
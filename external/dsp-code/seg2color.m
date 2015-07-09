% function to get the color map of object segmentation
function Im=seg2color(seg,objcolormap)

seg = seg +1;
[height,width]=size(seg);
Im=zeros(height*width,3);
nClasses=size(objcolormap,1);

for i=1:nClasses
    index=find(seg==i);
    if ~isempty(index)
        Im(index,:)=repmat(objcolormap(i,:),[length(index),1]);
    end
end

Im=reshape(Im,[height width 3]);
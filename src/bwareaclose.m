function bw2 = bwareaclose(bw,p,conn) 

CC = bwconncomp(bw,conn);
area = cellfun(@numel, CC.PixelIdxList);

idxToKeep = CC.PixelIdxList(area <= p);
idxToKeep = vertcat(idxToKeep{:});

bw2 = false(size(bw));
bw2(idxToKeep) = true;

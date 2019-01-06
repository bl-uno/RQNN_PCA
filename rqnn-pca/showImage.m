function showImage(W,k,sz,repet,maxrepet) 


for i=1:k
    subplot(maxrepet,k,(repet-1)*k+i)
    W1 =reshape(W(:,i),sz(1),sz(2));
    imshow(W1,[min(min(W1)) max(max(W1))])
end

function reference = find_ref(sig,fsamp)

for i = 1 : size(sig,1)
   e=1;
   for j = 1 : size(sig,1)
      if j~=i,
         [dist(e), scale, b]=sig_dis(sig(i,:),sig(j,:),fsamp,'dist_scale');
         e=e+1;
      end;
   end;
   aver_dist(i) = mean(dist);
   clear dist
end;

[a reference]=min(aver_dist);


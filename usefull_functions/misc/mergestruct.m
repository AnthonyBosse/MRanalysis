function [c] = mergestruct(a,b);
%  [am,bm] = mergestruct(a,b);
%
% merge 2 structures that are dissimilar with empty fields in both to make them match
% testor@locean-ipsl.upmc.fr, 2011/06/05
%
am = a;
bm = b;

     aaa = fieldnames(a);
     bbb = fieldnames(b);

     for jj = 1:length(aaa),
       if any(strcmp(aaa{jj},bbb))
       else
         for jk = 1:length(b)
         eval(['bm(jk).' aaa{jj} '= [];'])
         end
       end
     end

     for jj = 1:length(bbb),
       if any(strcmp(bbb{jj},aaa))
       else
         for jk = 1:length(a)
         eval(['am(jk).' bbb{jj} '= [];'])
         end
       end
     end

     c = [am bm];

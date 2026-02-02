seismo=load('./earthquakes.mat'); 
seismo=seismo.seismo;
ns=length(seismo.date);

todel=[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26];
%         for ja=length(evlo):-1:1
%             if(ismember(ja,todel))
%             evlo(ja)=[];
%             evla(ja)=[];
%             evdp(ja)=[];
%             end
%         end

for i=ns:-1:1
    if(ismember(i,todel))
        seismo.date(i)=[];
        seismo.time(i)=[];
        seismo.lat(i)=[];
        seismo.lon(i)=[];
        seismo.depth(i)=[];
        seismo.m(i)=[];
        seismo.distance(i)=[];
        seismo.az(i)=[];
    end
end

 namef=strcat('earthquakes1.mat');
 save(namef,'seismo');
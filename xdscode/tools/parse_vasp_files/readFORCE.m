function FC = readFORCE(file)
%%
% fp=fopen(file,'r'); 
% 
% nats=fscanf(fp,'%d',[2,1])
% FC = zeros(3,3,nats(1),nats(1));

% aa=fscanf(fp,'%d %d\n',[2,1])
% while(aa)
%     FC(:,:,aa(2),aa(1))=fscanf(fp,'%lf',[3,3]);
%     aa=fscanf(fp,'%d',[2,1])
% 
% end
% fclose(fp);

%%
fp=fopen(file,'r'); % open file descriptor
nats=fscanf(fp,'%d',[2,1]);
FC = zeros(3,3,nats(1),nats(1));
[a,m]=fscanf(fp,'%d',[2,1]);
while(m)
    FC(:,:,a(2),a(1))=fscanf(fp,'%lf',[3,3]);
    [a,m]=fscanf(fp,'%d',[2,1]);
end
fclose(fp);

end

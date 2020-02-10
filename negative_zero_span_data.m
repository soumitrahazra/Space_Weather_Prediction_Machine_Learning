clear all; clc;

afid=fopen('final_noflare_negative_harp_data3.csv');
T= table;

for i=0:1
line1=fgets(afid);
end

no=0;
  
for i=0:3766
  if line1 ~=-1
      if i < 28
       no=no+1;
   harp_no(no)= str2double(line1(1:2));
   line1=fgets(afid);
   
      else  if i < 506
   no=no+1;
   harp_no(no)= str2double(line1(1:3));
   line1=fgets(afid);
      else
     no=no+1;
   harp_no(no)= str2double(line1(1:5));
   line1=fgets(afid);     
          end  
      end
 end
end
fclose(afid)

for j=1:length(harp_no)
 j
nas=harp_no(j)
fname = sprintf('%d.txt', nas);
input_data= '/home/soumitra/Solar_Flare_Project_Machine_Learning/SHARPS';
sfid= fopen(fullfile(input_data,fname));
year=0;
mm=0;
dd=0;
hh=0;
mn=0;
ss=0;
times=0;

for i=0:1
line2=fgets(sfid);
end

nos=0;
while (line2 ~=-1)
nos=nos+1;
year(nos)=str2double(line2(1:4));
mm(nos)=str2double(line2(6:7));
dd(nos)=str2double(line2(9:10));
hh(nos)=str2double(line2(12:13));
mn(nos)=str2double(line2(15:16));
ss(nos)=str2double(line2(18:19));
Dates= [year(nos), mm(nos), dd(nos), hh(nos), mn(nos), ss(nos)];
final_dates(nos)= datetime(Dates);
times(nos)=str2double(final_dates(nos));
line2=fgets(sfid);
end
fclose(sfid);

js=fix(nos/2)
t_flare=final_dates(js);
tmax=t_flare-hours(24);
tmin=t_flare-hours(25);


so=0;
sa=[];
nop= length(final_dates);
for i=1:nos
 if (final_dates(i) >= tmin) && (final_dates(i) <= tmax)
  so=so+1;
  sa(so)=i;
 end
 end

n1=min(sa);
n2=max(sa);

%dlmread was used to read .txt file with one row and one column offset
%one can use csvread command if the file is in csv format.
aa=dlmread(fullfile(input_data,fname),'',1,1);
fa=size(aa);
ny=fa(2);

asf= aa(max(sa), 1:ny);
ast= final_dates(max(sa));
asl=(0.0*max(sa))/max(sa);
T_deg = table(ast', asf(:,1), asf(:,2), asf(:,3), asf(:,4), ...
asf(:,5), asf(:,6), asf(:,7), asf(:,8), asf(:,9), asf(:,10), ...
asf(:,11), asf(:,12), asf(:,13), asf(:,14), asf(:,15), asf(:,16), asl');
T = [T; T_deg];

end

T.Properties.VariableNames = {'Date_Obs', 'USFLUX','MEANGAM','MEANGBT','MEANGBZ','MEANGBH', ...
 'MEANJZD','TOTUSJZ','MEANALP','MEANJZH','TOTUSJH','ABSNJZH','SAVNCPP', ...
 'MEANPOT','TOTPOT','MEANSHR','SHRGT45', 'Label'};

filename = sprintf('negative_24_hour_back_noflare.csv');
%Writing Csv file into another folder
mypath='/home/soumitra/Flare_Paper_Final_Result/Another_Data_Set_loop24span0/Final_Combine';
%Lets write the table in a csv file
writetable(T,fullfile(mypath,filename),'Delimiter',',');



clear all; clc;

afid=fopen('final_flare_positive_harp_data3.csv');
T=table;

for i=0:1
line1=fgets(afid);
end

no=0;
  for i=0:690
  if line1 ~=-1
   if i < 91
   no=no+1;
   ar_no(no)=str2double(line1(6:10));
   harp_no(no)= str2double(line1(11:14));
   yr(no)=str2double(line1(17:21));
   mma(no)=str2double(line1(23:24));
   day(no)=str2double(line1(26:27));
    hr(no)=str2double(line1(29:30));
   mins(no)= str2double(line1(32:33));
   ssa(no)=0.0;
   Date= [yr(no), mma(no), day(no), hr(no), mins(no), ssa(no)];
   final_date(no)= datetime(Date);
   line1=fgets(afid);
      else
     no=no+1;
   ar_no(no)=str2double(line1(6:10));
   harp_no(no)= str2double(line1(11:15));
   yr(no)=str2double(line1(18:22));
   mma(no)=str2double(line1(24:25));
   day(no)=str2double(line1(27:28));
    hr(no)=str2double(line1(30:31));
   mins(no)= str2double(line1(33:34));
   ssa(no)=0.0;
   Date= [yr(no), mma(no), day(no), hr(no), mins(no), ssa(no)];
   final_date(no)= datetime(Date);
   line1=fgets(afid);     
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

t_flare=final_date(j);
tmax=t_flare-hours(24);
tmin=t_flare-hours(25);

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

so=0;
sa=[];
nop= length(final_dates);
for i=1:nos
 if (final_dates(i) >= tmin) && (final_dates(i) <= tmax)
  so=so+1;
  sa(so)=i;
 end
end
 
%dlmread was used to read .txt file with one row and one column offset
%one can use csvread command if the file is in csv format.
aa=dlmread(fullfile(input_data,fname),'',1,1);
fa=size(aa);
ny=fa(2);
asf= aa(max(sa), 1:ny);
ast= final_dates(max(sa));
asl=max(sa)/max(sa);
T_deg = table(ast', asf(:,1), asf(:,2), asf(:,3), asf(:,4), ...
asf(:,5), asf(:,6), asf(:,7), asf(:,8), asf(:,9), asf(:,10), ...
asf(:,11), asf(:,12), asf(:,13), asf(:,14), asf(:,15), asf(:,16), asl');
T = [T; T_deg];

end
T.Properties.VariableNames = {'Date_Obs', 'USFLUX','MEANGAM','MEANGBT','MEANGBZ','MEANGBH', ...
 'MEANJZD','TOTUSJZ','MEANALP','MEANJZH','TOTUSJH','ABSNJZH','SAVNCPP', ...
 'MEANPOT','TOTPOT','MEANSHR','SHRGT45', 'Label'};

filename = sprintf('positive_flare_24_hour_back.csv');
%Writing Csv file into another folder
mypath='/home/soumitra/Flare_Paper_Final_Result/Another_Data_Set_loop24span0/Final_Combine';
%Lets write the table in a csv file
writetable(T,fullfile(mypath,filename),'Delimiter',',');




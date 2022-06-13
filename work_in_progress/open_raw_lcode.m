O = OsirisLoader();

% openID = fopen(['r2l_2x','/',obj.field,obj.direction,lcode_dump_char,'.swp'],'r');
lcode_dump_char = '00005';
raw_name = 'tb'; 
openID = fopen(['r2l_2x','/',raw_name,lcode_dump_char,'.swp'],'r');

data_temp = fread(openID,[8,inf],'double')';

% z,r,pz,pr,pa,q,w,id

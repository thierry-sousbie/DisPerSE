CHECK_ALL_FILES = 0;
LOAD_CALIBOBJID = 0;
LOAD_DATA = 0;

photo_rows = [30,31,32,33];
photo_rows_name = ["Q","Qerr","U","Uerr"];
photo_filters = ["u","g","r","i","z"];

spec_rows = [4,5,6];
spec_rows_name = ["ra","dec","z"];


addr="http://sdss.physics.nyu.edu/vagc-dr7/vagc2/sdss/parameters/";
base="calibObj-";
fitsdir="./photo/";

LIGHTSPEED = 299792.458;


if (LOAD_CALIBOBJID) {
  u=fits_read("object_sdss_imaging.fits",h,hdu=2);
  calibObj_ID = (*u(33))(int(id(1,))+1);
  u=[];
  save,createb("calibObjId.pdb"),calibObj_ID;
 } else restore,openb("calibObjId.pdb");
 

lss=read_ascii("lss.dr72bright0.dat");lss(6,)/=LIGHTSPEED;
id=read_ascii("id.dr72bright0.dat");
run=int(id(2,));
camcol=int(id(4,));

uniqueID=unique(run*10+camcol);
if (LOAD_DATA) datap=array(pointer,max(uniqueID));
r=uniqueID/10;
c=uniqueID-10*r;
N=0;

for (i=1;i<=numberof(uniqueID);i++) {
  fname = swrite(format="%s%6.6d-%d.fits",base,r(i),c(i));
  write,format="%d/%d: %s \n",i,numberof(uniqueID),fname;
  if ((CHECK_ALL_FILES)||(!fileExist(fitsdir+fname))) {
    write,format="CHEKING (%d)\n",++N;
    exec,"wget "+addr+fname+" -N -P "+fitsdir;
  }
  
  if (LOAD_DATA) {
    u=fits_read(fitsdir+fname,h,hdu=2);
    res=[];
    for (j=1;j<=numberof(photo_rows);j++) grow,res,[*u(photo_rows(j))];
    datap(uniqueID(i))=&res;
  }
 }

if (LOAD_DATA) save,createb("photo_data.pdb"),datap; else restore,openb("photo_data.pdb");


//r=[];c=[];
old_data_id = -1;

result=array(double,[2,numberof(spec_rows)+numberof(photo_filters)*numberof(photo_rows),dimsof(lss)(0)]);
result(1:numberof(spec_rows),)=lss(spec_rows,);

for (i=1;i<=numberof(run);i++) {
  data_id = run(i)*10+camcol(i);
  //error,"stop";
  //if (old_data_id != data_id)
  data = *datap(data_id);
  //old_data_id = data_id;
  
  result(numberof(spec_rows)+1:,i)=transpose(data(calibObj_ID(i)+1,,))(*);
 }
  
header = spec_rows_name(1);
for (i=2;i<=numberof(spec_rows);i++) header=swrite(format="%s %s",header,spec_rows_name(i));

for (j=1;j<=numberof(photo_filters);j++)
  for (i=1;i<=numberof(photo_rows);i++)
    header=swrite(format="%s %s",header,photo_rows_name(i)+"_"+photo_filters(j));

smwrite("catalog.dat",transpose(result),head=header);


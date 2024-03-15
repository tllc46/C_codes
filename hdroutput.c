#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_f(char *field,float val)
{
	if(val==-12345)
	{
		printf("%s = undef\n",field);
	}
	else
	{
		printf("%s = %f\n",field,val);
	}
}

void print_n(char *field,int val)
{
	if(val==-12345)
	{
		printf("%s = undef\n",field);
	}
	else
	{
		printf("%s = %d\n",field,val);
	}
}

void print_k(char *field,char *val)
{
	char undef[7]="-12345";
	if(!strncmp(val,undef,6))
	{
		printf("%s = undef\n",field);
	}
	else
	{	
		char dum[9]="   kevnm";
		if(!strncmp(field,dum,8))
		{
			char str[17];
			strncpy(str,val,16);
			str[16]='\0';
			printf("%s = %s\n",field,str);
		}
		else
		{
			char str[9];
			strncpy(str,val,8);
			str[8]='\0';
			printf("%s = %s\n",field,str);
		}
	}
}

void print_f64(char *field,double val)
{
	if(val==-12345)
	{
		printf("%s = undef\n",field);
	}
	else
	{
		printf("%s = %f\n",field,val);
	}
}	

int main(int argc,char **argv)
{
	FILE *fp;
	float f[70];
	int n[40];
	char k[192];
	int e;
	size_t b;

	fp=fopen(argv[1],"rb");
	b=fread(f,sizeof(float),70,fp);
	b=fread(n,sizeof(int),40,fp);
	b=fread(k,sizeof(char),192,fp);

	printf("<required fields>\n");
	print_n("    npts",n[9]);
	print_n("   nvhdr",n[6]);
	print_f("       b",f[5]);
	print_f("       e",f[6]);
	print_n("  iftype",n[15]);
	print_n("   leven",n[35]);
	print_f("   delta",f[0]);

	printf("\n");

	printf("<time fields>\n");
	print_f("  odelta",f[4]);
	print_n("    idep",n[16]);
	print_f("  depmin",f[1]);
	print_f("  depmax",f[2]);
	print_f("  depmen",f[56]);
	print_n("  nzyear",n[0]);
	print_n("  nzjday",n[1]);
	print_n("  nzhour",n[2]);
	print_n("   nzmin",n[3]);
	print_n("   nzsec",n[4]);
	print_n("  nzmsec",n[5]);
	print_n("  iztype",n[17]);
	print_f("       o",f[7]);
	print_k("      ko",k+32);
	print_n("  nsnpts",n[10]);
	print_f("      sb",f[54]);
	print_f("  sdelta",f[55]);

	printf("\n");

	printf("<phase picks>\n");
	print_f("       a",f[8]);
	print_k("      ka",k+40);
	print_f("       f",f[20]);
	print_k("      kf",k+128);
	print_f("      t0",f[10]);
	print_f("      t1",f[11]);
	print_f("      t2",f[12]);
	print_f("      t3",f[13]);
	print_f("      t4",f[14]);
	print_f("      t5",f[15]);
	print_f("      t6",f[16]);
	print_f("      t7",f[17]);
	print_f("      t8",f[18]);
	print_f("      t9",f[19]);
	print_k("     kt0",k+48);
	print_k("     kt1",k+56);
	print_k("     kt2",k+64);
	print_k("     kt3",k+72);
	print_k("     kt4",k+80);
	print_k("     kt5",k+88);
	print_k("     kt6",k+96);
	print_k("     kt7",k+104);
	print_k("     kt8",k+112);
	print_k("     kt9",k+120);

	printf("\n");

	printf("<instrument fields>\n");
	print_k("   kinst",k+184);
	print_n("   iinst",n[19]);
	print_f("   resp0",f[21]);
	print_f("   resp1",f[22]);
	print_f("   resp2",f[23]);
	print_f("   resp3",f[24]);
	print_f("   resp4",f[25]);
	print_f("   resp5",f[26]);
	print_f("   resp6",f[27]);
	print_f("   resp7",f[28]);
	print_f("   resp8",f[29]);
	print_f("   resp9",f[30]);

	printf("\n");

	printf("<station fields>\n");
	print_k("  knetwk",k+168);
	print_k("   kstnm",k);
	print_n("  istreg",n[20]);
	print_f("    stla",f[31]);
	print_f("    stlo",f[32]);
	print_f("    stel",f[33]);
	print_f("    stdp",f[34]);
	print_f("   cmpaz",f[57]);
	print_f("  cmpinc",f[58]);
	print_k("  kcmpnm",k+160);
	print_n("  lpspol",n[36]);

	printf("\n");

	printf("<event fields>\n");
	print_k("   kevnm",k+8);
	print_n("  ievreg",n[21]);
	print_f("    evla",f[35]);
	print_f("    evlo",f[36]);
	print_f("    evel",f[37]);
	print_f("    evdp",f[38]);
	print_f("     mag",f[39]);
	print_n(" imagtyp",n[25]);
	print_n(" imagsrc",n[26]);
	print_n("  ievtyp",n[22]);
	print_n("   nevid",n[8]);
	print_n("   norid",n[7]);
	print_n("   nwfid",n[11]);
	print_k("   khole",k+24);
	print_f("    dist",f[50]);
	print_f("      az",f[51]);
	print_f("     baz",f[52]);
	print_f("   gcarc",f[53]);
	print_n("   ibody",n[27]);

	printf("\n");

	printf("<miscellaneous fields>\n");
	print_n("  lcalda",n[38]);
	print_n("   iqual",n[23]);
	print_n("  isynth",n[24]);
	print_k("  kdatrd",k+176);
	print_f("   user0",f[40]);
	print_f("   user1",f[41]);
	print_f("   user2",f[42]);
	print_f("   user3",f[43]);
	print_f("   user4",f[44]);
	print_f("   user5",f[45]);
	print_f("   user6",f[46]);
	print_f("   user7",f[47]);
	print_f("   user8",f[48]);
	print_f("   user9",f[49]);
	print_k("  kuser0",k+136);
	print_k("  kuser1",k+144);
	print_k("  kuser2",k+152);
	print_n("  lovrok",n[37]);
	print_n("  nxsize",n[12]);
	print_n("  nysize",n[13]);
	print_f("xminimum",f[59]);
	print_f("xmaximum",f[60]);
	print_f("yminimum",f[61]);
	print_f("ymaximum",f[62]);

	printf("\n");

	printf("<unused fields>\n");
	print_f("   scale",f[3]);
	print_f("     fmt",f[9]);
	print_f("   adjtm",f[63]);
	print_f("  fhdr65",f[64]);
	print_f("  fhdr66",f[65]);
	print_f("  fhdr67",f[66]);
	print_f("  fhdr68",f[67]);
	print_f("  fhdr69",f[68]);
	print_f("  fhdr70",f[69]);
	print_n("  nhdr15",n[14]);
	print_n("   ihdr4",n[18]);
	print_n("  ihdr14",n[28]);
	print_n("  ihdr15",n[29]);
	print_n("  ihdr16",n[30]);
	print_n("  ihdr17",n[31]);
	print_n("  ihdr18",n[32]);
	print_n("  ihdr19",n[33]);
	print_n("  ihdr20",n[34]);
	print_n("   lhdr5",n[39]);

	if(n[6]==7)
	{
		double f64[22];
		e=fseek(fp,n[9]*sizeof(float),SEEK_CUR);
		b=fread(f64,sizeof(double),22,fp);
		printf("\n");
		
		printf("<footer fields>\n");
		print_f64("     b",f64[1]);
		print_f64("     e",f64[2]);
		print_f64(" delta",f64[0]);
		print_f64("     o",f64[3]);
		print_f64("    sb",f64[20]);
		print_f64("sdelta",f64[21]);
		print_f64("     a",f64[4]);
		print_f64("     f",f64[15]);
		print_f64("    t0",f64[5]);
		print_f64("    t1",f64[6]);
		print_f64("    t2",f64[7]);
		print_f64("    t3",f64[8]);
		print_f64("    t4",f64[9]);
		print_f64("    t5",f64[10]);
		print_f64("    t6",f64[11]);
		print_f64("    t7",f64[12]);
		print_f64("    t8",f64[13]);
		print_f64("    t9",f64[14]);
		print_f64("  stla",f64[19]);
		print_f64("  stlo",f64[18]);
		print_f64("  evla",f64[17]);
		print_f64("  evlo",f64[16]);
	}	

	e=fclose(fp);
	return 0;
}

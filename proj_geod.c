/*
gcc proj_geod.c -L/usr/lib/x86_64-linux-gnu -lproj
*/

#include <stdio.h>
#include <proj.h>

int main(int argc,char **argv)
{
        PJ *obj;
        PJ_PROJ_INFO info;
        PJ_COORD coord1;
        PJ_COORD coord2;
        PJ_COORD distaz;

        obj=proj_create(PJ_DEFAULT_CTX,"+proj=lonlat +ellps=WGS84");

        info=proj_pj_info(obj);
        printf("id=%s\n",info.id);
        printf("description=%s\n",info.description);
        printf("definition=%s\n",info.definition);
        printf("has_inverse=%d\n",info.has_inverse);
        printf("accuracy=%f\n",info.accuracy);
        printf("\n");

        coord1.lp.lam=proj_torad(0);
        coord1.lp.phi=proj_torad(0);
        coord2.lp.lam=proj_torad(10);
        coord2.lp.phi=proj_torad(10);

        distaz=proj_geod(obj,coord1,coord2);

        printf("distance=%f\n",distaz.geod.s);
        printf("azimuth=%f\n",distaz.geod.a1);
        printf("back azimuth=%f\n",distaz.geod.a2);

        proj_destroy(obj);

        return 0;
}

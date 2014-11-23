#include <stdio.h>
#include <stdlib.h>
#include "coordtrans.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s from to lng lat\n", argv[0]);
        return 1;
    }

    double lng = atof(argv[3]);
    double lat = atof(argv[4]);
    double newlng = 0;
    double newlat = 0;
    int ret = coordtrans(argv[1], argv[2], lng, lat, newlng, newlat);

    if (ret == 0) {
        printf("%lf\t%lf\t%lf,%lf\n", newlng, newlat, newlng, newlat);
        return 0;
    } else {
        printf("failed to convert from %s to %s\n", argv[1], argv[2]);
        return 1;
    }
}

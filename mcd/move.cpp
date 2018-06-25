#include "stdafx.h"

struct wallcre {
	int type;
	double x1,y1,z1,x2,y2,z2;
	double depth, area;
};

struct wallbre {
	double x1, y1, z1, x2, y2, z2;
};

wallcre wallcarray[120];
wallbre wallbarray[40];

int wallc(movestate *pe, int(*block)[256][256]) {
	int i, j, k, l;
	l = 0;
	int wallbn, wallcn;
	wallbn = 0;
	for (k = (int)(pe->pos.x - r - rm); k <= (int)(pe->pos.x + r + rm); k++) {
		for (i = (int)(pe->pos.z - r - rm); i <= (int)(pe->pos.z + r + rm); i++) {
			for (j = (int)(pe->pos.y + bot - rm); j <= (int)(pe->pos.y + top + rm); j++) {
				if ((0 <= i)&(i <= 255)&(0 <= j)&(j <= 255)&(0 <= k)&(k <= 255)) {
					if (!tr[block[i][j][k] & 0x00ff]) {
						wallbarray[wallbn].x1 = k;
						wallbarray[wallbn].x2 = k + 1;
						wallbarray[wallbn].y1 = j;
						wallbarray[wallbn].y2 = j + 1;
						wallbarray[wallbn].z1 = i;
						wallbarray[wallbn].z2 = i + 1;
						wallbn++;
					}
				}
			}
		}
	}
	for (i = 0; i < 4; i++) {
		wallcn = 0;
		k = 0;
		for (j = 0; j < wallbn; j++) {
			if (pe->pos.x >(wallbarray[j].x1 - r) && pe->pos.x < wallbarray[j].x1 && (pe->pos.y + bot) < wallbarray[j].y2 && (pe->pos.y + top) > wallbarray[j].y1 && pe->pos.z<(wallbarray[j].z2 + r) && pe->pos.z>(wallbarray[j].z1 - r)) {
				wallcarray[wallcn].type = 0;
				wallcarray[wallcn].depth = wallbarray[j].x1 - r;
				wallcarray[wallcn].area = (min(pe->pos.z + r, wallbarray[j].z2) - max(pe->pos.z - r, wallbarray[j].z1))*(min(pe->pos.y + top, wallbarray[j].y2) - max(pe->pos.y + bot, wallbarray[j].y1));
				wallcn++;
			}
			if (pe->pos.x < (wallbarray[j].x2 + r) && pe->pos.x > wallbarray[j].x2 && (pe->pos.y + bot) < wallbarray[j].y2 && (pe->pos.y + top) > wallbarray[j].y1 && pe->pos.z<(wallbarray[j].z2 + r) && pe->pos.z>(wallbarray[j].z1 - r)) {
				wallcarray[wallcn].type = 1;
				wallcarray[wallcn].depth = wallbarray[j].x2 + r;
				wallcarray[wallcn].area = (min(pe->pos.z + r, wallbarray[j].z2) - max(pe->pos.z - r, wallbarray[j].z1))*(min(pe->pos.y + top, wallbarray[j].y2) - max(pe->pos.y + bot, wallbarray[j].y1));
				wallcn++;
			}
			if (pe->pos.z > (wallbarray[j].z1 - r) && pe->pos.z < wallbarray[j].z1 && (pe->pos.y + bot) < wallbarray[j].y2 && (pe->pos.y + top) > wallbarray[j].y1 && pe->pos.x<(wallbarray[j].x2 + r) && pe->pos.x>(wallbarray[j].x1 - r)) {
				wallcarray[wallcn].type = 4;
				wallcarray[wallcn].depth = wallbarray[j].z1 - r;
				wallcarray[wallcn].area = (min(pe->pos.x + r, wallbarray[j].x2) - max(pe->pos.x - r, wallbarray[j].x1))*(min(pe->pos.y + top, wallbarray[j].y2) - max(pe->pos.y + bot, wallbarray[j].y1));
				wallcn++;
			}
			if (pe->pos.z < (wallbarray[j].z2 + r) && pe->pos.z > wallbarray[j].z2 && (pe->pos.y + bot) < wallbarray[j].y2 && (pe->pos.y + top) > wallbarray[j].y1 && pe->pos.x<(wallbarray[j].x2 + r) && pe->pos.x>(wallbarray[j].x1 - r)) {
				wallcarray[wallcn].type = 5;
				wallcarray[wallcn].depth = wallbarray[j].z2 + r;
				wallcarray[wallcn].area = (min(pe->pos.x + r, wallbarray[j].x2) - max(pe->pos.x - r, wallbarray[j].x1))*(min(pe->pos.y + top, wallbarray[j].y2) - max(pe->pos.y + bot, wallbarray[j].y1));
				wallcn++;
			}
			if ((pe->pos.y + top) > wallbarray[j].y1 && pe->pos.y < wallbarray[j].y1 && pe->pos.x<(wallbarray[j].x2 + r) && pe->pos.x>(wallbarray[j].x1 - r) && pe->pos.z<(wallbarray[j].z2 + r) && pe->pos.z>(wallbarray[j].z1 - r)) {
				wallcarray[wallcn].type = 2;
				wallcarray[wallcn].depth = wallbarray[j].y1 - top;
				wallcarray[wallcn].area = (min(pe->pos.x + r, wallbarray[j].x2) - max(pe->pos.x - r, wallbarray[j].x1))*(min(pe->pos.z + r, wallbarray[j].z2) - max(pe->pos.z - r, wallbarray[j].z1));
				wallcn++;
			}
			if ((pe->pos.y + bot) < wallbarray[j].y2 && pe->pos.y > wallbarray[j].y2 && pe->pos.x<(wallbarray[j].x2 + r) && pe->pos.x>(wallbarray[j].x1 - r) && pe->pos.z<(wallbarray[j].z2 + r) && pe->pos.z>(wallbarray[j].z1 - r)) {
				wallcarray[wallcn].type = 3;
				wallcarray[wallcn].depth = wallbarray[j].y2 - bot;
				wallcarray[wallcn].x1 = wallbarray[j].x1;
				wallcarray[wallcn].z1 = wallbarray[j].z1;
				wallcarray[wallcn].y1 = wallbarray[j].y1;
				wallcarray[wallcn].area = (min(pe->pos.x + r, wallbarray[j].x2) - max(pe->pos.x - r, wallbarray[j].x1))*(min(pe->pos.z + r, wallbarray[j].z2) - max(pe->pos.z - r, wallbarray[j].z1));
				wallcn++;
			}
		}
		if (wallcn == 0)
			break;
		for (j = 0; j < wallcn; j++) {
			if (wallcarray[j].area > wallcarray[k].area) {
				k = j;
			}
		}
		switch (wallcarray[k].type) {
		case 0:
		case 1:
			pe->pos.x = wallcarray[k].depth;
			pe->v.x = 0;
			break;
		case 3:
			pe->flying = 0;
			l = 1;
			pe->onlandstate.x1 = wallcarray[k].x1;
			pe->onlandstate.z1 = wallcarray[k].z1;
			pe->onlandstate.y = wallcarray[k].y1;
		case 2:
			pe->pos.y = wallcarray[k].depth;
			pe->v.y = 0;
			break;
		case 5:
		case 4:
			pe->pos.z = wallcarray[k].depth;
			pe->v.z = 0;
			break;
		}
	}
	if (pe->onlandstate.isonland && (l == 0) && (pe->up == -1) && (block[(int)(pe->onlandstate.z1 + 0.5)][(int)(pe->onlandstate.y + 0.5)][(int)(pe->onlandstate.x1 + 0.5)])) {
		l = 1;
		pe->pos.x = pe->pos.x > (pe->onlandstate.x1 - r + 0.001) ? (pe->pos.x) : ((pe->v.x = 0), pe->onlandstate.x1 - r + 0.001);
		pe->pos.z = pe->pos.z > (pe->onlandstate.z1 - r + 0.001) ? (pe->pos.z) : ((pe->v.z = 0), pe->onlandstate.z1 - r + 0.001);
		pe->pos.x = pe->pos.x < (pe->onlandstate.x1 + 1 + r - 0.001) ? (pe->pos.x) : ((pe->v.x = -0), pe->onlandstate.x1 + 1 + r - 0.001);
		pe->pos.z = pe->pos.z < (pe->onlandstate.z1 + 1 + r - 0.001) ? (pe->pos.z) : ((pe->v.z = -0), pe->onlandstate.z1 + 1 + r - 0.001);
		pe->pos.y = pe->onlandstate.y - bot + 1;
		pe->v.y = 0;
	}
	pe->onlandstate.isonland = l;
	return 0;
}

int move(movestate *pe, unsigned int dt, int(*block)[256][256]) {
	double a0x = 0.0, a0y = 0.0, a0z = 0.0;
	double a, v;
	int i, j;
	if (dt > 1000)
		return 1;
	for (i = 0; i < dt; i++) {
		a = sqrt(pe->forward*pe->forward + pe->left * pe->left);
		if (pe->flying) {
			a0y = pe->up * 30.0 - pe->v.y * 5.0;
		}
		else {
			a0y = -20;
		}
		if (pe->forward != 0 || pe->left != 0) {
			a0x = (cos(pe->ang1)*pe->forward + sin(pe->ang1) * pe->left)*40.0 / a;
			a0z = (sin(pe->ang1)*pe->forward - cos(pe->ang1) * pe->left)*40.0 / a;
		}
		pe->v.x += a0x / 1000;
		pe->v.y += a0y / 1000;
		pe->v.z += a0z / 1000;
		v = sqrt(pe->v.x*pe->v.x + pe->v.z * pe->v.z);
		if (pe->flying) {
			if (v > 0.005) {
				pe->v.x -= (3 * pe->v.x + 4 * pe->v.x / v) / 1000;
				pe->v.z -= (3 * pe->v.z + 4 * pe->v.z / v) / 1000;
			}
			else {
				pe->v.z = 0;
				pe->v.x = 0;
			}
		}
		else {
			if (v > 0.021) {
				pe->v.x -= (5 * pe->v.x + 15 * pe->v.x / v) / 1000;
				pe->v.z -= (5 * pe->v.z + 15 * pe->v.z / v) / 1000;
			}
			else {
				pe->v.z = 0;
				pe->v.x = 0;
			}
		}
		pe->pos.x += pe->v.x / 1000;
		pe->pos.y += pe->v.y / 1000;
		pe->pos.z += pe->v.z / 1000;
		wallc(pe, block);
	}
	return 0;
}


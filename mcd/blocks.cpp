#include "stdafx.h"

int upbs = 4;

int *bneedup1;
int *bneedup2;
int maxupbn;

int updio(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate);
int uprt(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate);
int upitb(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate);
int uprl(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate);
int glev(int x, int y, int z, int f, int(*block)[256][256]);
int glevdio(int x, int y, int z, int f, int(*block)[256][256]);
void addtolist(int x, int y, int z);
void removefromlist(int x, int y, int z);

void addtolist(int x, int y, int z) {
	int i;
	for (i = 0; i < maxupbn; i++) {
		if (bneedup1[i] == x + (y << 8) + (z << 16) + 0x01000000) {
			break;
		}
	}
	if (i == maxupbn) {
		for (i = 0; i < maxupbn; i++) {
			if (bneedup1[i] == 0) {
				bneedup1[i] = x + (y << 8) + (z << 16) + 0x01000000;
				break;
			}
		}
		if (i == maxupbn) {
			int *t, j;
			t = bneedup1;
			bneedup1 = (int*)malloc(sizeof(int)*maxupbn * 2);
			for (j = 0; j < maxupbn; j++) {
				bneedup1[j] = t[j];
				bneedup1[j + maxupbn] = 0;
			}
			free(t);
			bneedup1[maxupbn] = x + (y << 8) + (z << 16) + 0x01000000;
			t = bneedup2;
			bneedup2 = (int*)malloc(sizeof(int)*maxupbn * 2);
			for (j = 0; j < maxupbn; j++) {
				bneedup2[j] = t[j];
				bneedup2[j + maxupbn] = 0;
			}
			free(t);
			maxupbn *= 2;
		}
	}
}

void removefromlist(int x, int y, int z) {
	int i;
	for (i = 0; i < maxupbn; i++) {
		if (bneedup1[i] == x + (y << 8) + (z << 16) + 0x01000000) {
			bneedup1[i] = 0;
			break;
		}
	}
}

void updateintick(int(*block)[256][256], elementtypes *regionsupdate) {
	int j;
	upbs = (upbs == 3) ? 4 : 3;
	for (j = 0; j < maxupbn; j++) {
		bneedup2[j] = bneedup1[j];
		bneedup1[j] = 0;
	}
	for (j = 0; j < maxupbn; j++) {
		if (bneedup2[j] & 0x01000000) {
			upb(bneedup2[j] & 0x00ff, (bneedup2[j] & 0xff00) >> 8, (bneedup2[j] & 0x00ff0000) >> 16, 3, block, regionsupdate);
		}
	}
}

int initblocks(int (**block)[256][256],double3 pos) {
	*block = (int(*)[256][256])malloc(sizeof(int) * 256 * 256 * 256);
	int i, j, k;
	for (i = 0; i<256; i++) {
		for (j = 0; j<256; j++) {
			for (k = 0; k<256; k++)
				(*block)[i][j][k] = 0;
		}
	}
	for (i = 1; i < 255; i++) {
		for (k = 1; k < 255; k++) {
			for (j = 1; j < 126; j++) {
				(*block)[i][j][k] = 1;
			}
			for (j = 126; j < 128; j++) {
				(*block)[i][j][k] = 1; 
			}
			(*block)[i][j][k] = 1;
		}
	}
	for (i = pos.z - r - 0.1; i < pos.z + 0.1 + r; i++) {
		for (k = pos.x - r - 0.1; k <pos.x + 0.1 + r; k++) {
			for (j = pos.y + bot - 0.1; j < pos.y + 0.1 + top; j++) {
				(*block)[i][j][k] = 0;
			}
		}
	}
	maxupbn = 64;
	bneedup1 = (int*)malloc(sizeof(int)*maxupbn);
	bneedup2 = (int*)malloc(sizeof(int)*maxupbn);
	for (i = 0; i < maxupbn; i++) {
		bneedup1[i] = 0;
		bneedup2[i] = 0;
	}
	return 0;
}

void removeblock(movestate *pe, int(*block)[256][256], elementtypes *regionsupdate) {
	int blockid, x, y, z;
	x = pe->point.blockx;
	y = pe->point.blocky;
	z = pe->point.blockz;
	blockid = block[z][pe->point.blocky][x];
	block[pe->point.blockz][pe->point.blocky][pe->point.blockx] = 0;
	regionsupdate[pe->point.blockx / 16 + pe->point.blockz / 16 * 16].opaque |= 1;
	if (!tr[blockid & 0xff]) {
		if (x % 16 == 0 && x > 0)
			regionsupdate[x / 16 + z / 16 * 16 - 1].opaque |= 1;
		if (z % 16 == 0 && z > 0)
			regionsupdate[x / 16 + z / 16 * 16 - 16].opaque |= 1;
		if (x % 16 == 15 && x < 255)
			regionsupdate[x / 16 + z / 16 * 16 + 1].opaque |= 1;
		if (z % 16 == 15 && z < 255)
			regionsupdate[x / 16 + z / 16 * 16 + 16].opaque |= 1;
	}
	if (x > 0) {
		upb(x - 1, y, z, 0, block, regionsupdate);
	}
	if (y > 0) {
		upb(x, y - 1, z, 0, block, regionsupdate);
	}
	if (z > 0) {
		upb(x, y, z - 1, 0, block, regionsupdate);
	}
	if (x < 255) {
		upb(x + 1, y, z, 0, block, regionsupdate);
	}
	if (y < 255) {
		upb(x, y + 1, z, 0, block, regionsupdate);
	}
	if (z < 255) {
		upb(x, y, z + 1, 0, block, regionsupdate);
	}
}

void rightbutton(movestate *pe, int(*block)[256][256], elementtypes *regionsupdate) {
	int x, y, z;
	x = pe->point.surx;
	y = pe->point.sury;
	z = pe->point.surz;
	if (pe->point.ispoint) {
		if ((block[pe->point.blockz][pe->point.blocky][pe->point.blockx] & 0x00ff) == 6 && pe->up != -1) {
			x = pe->point.blockx;
			y = pe->point.blocky;
			z = pe->point.blockz;
			block[z][y][x] = (block[z][y][x] & 0x03f3ff) + (((block[z][y][x] & 0x0c00) + 0x0400) & 0x0c00);
			regionsupdate[x / 16 + z / 16 * 16].opaque |= 1;
		}
		else if ((x >= 0) && (x <= 255) && (y >= 0) && (y <= 255) && (z >= 0) && (z <= 255) && !(pe->pos.x + r > x&&pe->pos.x - r<x + 1 && pe->pos.z + r>z&&pe->pos.z - r<z + 1 && pe->pos.y + top>y&&pe->pos.y + bot < y + 1 && tr[pe->inhand] == 0) && (block[z][y][x] & 0x00ff) == 0) {
			block[z][y][x] = pe->inhand;
			switch (pe->inhand) {
			case 5:
				block[z][y][x] |= 0x0800;
				if (pe->point.blocky == y - 1) {
					block[z][y][x] |= 0x0000;
				}
				else if (pe->point.blockx == x + 1) {
					block[z][y][x] |= 0x0100;
				}
				else if (pe->point.blockx == x - 1) {
					block[z][y][x] |= 0x0200;
				}
				else if (pe->point.blockz == z + 1) {
					block[z][y][x] |= 0x0300;
				}
				else if (pe->point.blockz == z - 1) {
					block[z][y][x] |= 0x0400;
				}
				else {
					block[z][y][x] = 0;
				}
				break;
			case 6: {
				if (pe->at.x >= pe->at.z && pe->at.x >= -pe->at.z) {
					block[z][y][x] |= 0x0000;
				}
				else if (-pe->at.x >= pe->at.z && -pe->at.x >= -pe->at.z) {
					block[z][y][x] |= 0x0100;
				}
				else if (pe->at.z >= pe->at.x && pe->at.z >= -pe->at.x) {
					block[z][y][x] |= 0x0200;
				}
				else if (-pe->at.z >= pe->at.x && -pe->at.z >= -pe->at.x) {
					block[z][y][x] |= 0x0300;
				}
				break;
			}
			}
			regionsupdate[x / 16 + z / 16 * 16].opaque |= 1;
			if (!tr[pe->inhand]) {
				if (x % 16 == 0 && x > 0)
					regionsupdate[x / 16 + z / 16 * 16 - 1].opaque |= 1;
				if (z % 16 == 0 && z > 0)
					regionsupdate[x / 16 + z / 16 * 16 - 16].opaque |= 1;
				if (x % 16 == 15 && x < 255)
					regionsupdate[x / 16 + z / 16 * 16 + 1].opaque |= 1;
				if (z % 16 == 15 && z < 255)
					regionsupdate[x / 16 + z / 16 * 16 + 16].opaque |= 1;
			}
		}
		upb(x, y, z, 0, block, regionsupdate);
		if (x > 0) {
			upb(x - 1, y, z, 0, block, regionsupdate);
		}
		if (y > 0) {
			upb(x, y - 1, z, 0, block, regionsupdate);
		}
		if (z > 0) {
			upb(x, y, z - 1, 0, block, regionsupdate);
		}
		if (x < 255) {
			upb(x + 1, y, z, 0, block, regionsupdate);
		}
		if (y < 255) {
			upb(x, y + 1, z, 0, block, regionsupdate);
		}
		if (z < 255) {
			upb(x, y, z + 1, 0, block, regionsupdate);
		}
	}
}

void upb(int x, int y, int z, int type, int(*block)[256][256], elementtypes *regionsupdate) {
	int ib;
	ib = block[z][y][x];
	switch (ib & 0x00ff) {
	case 0:
		return;
	case 4:
		uprl(x, y, z, type, block, regionsupdate);
		break;
	case 5:
		uprt(x, y, z, type, block, regionsupdate);
		break;
	case 6:
		updio(x, y, z, type, block, regionsupdate);
		break;
	default:
		upitb(x, y, z, type, block, regionsupdate);
	}
}

int updio(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate) {
	int ib, i, a, b, c, x1, z1, f;
	ib = block[z][y][x];
	if (s == 2) {
		block[z][y][x] = 0;
	}
	else if (y == 0) {
		block[z][y][x] = 0;
	}
	else if (tr[block[z][y - 1][x] & 0x00ff]) {
		block[z][y][x] = 0;
	}
	else {
		if (block[z][y][x] & ((upbs == 3) ? 0x020000 : 0x008000)) {
			if (true) {
				block[z][y][x] = block[z][y][x] & 0x017fff; //clear change in 3or4
				a = (block[z][y][x] & 0x003000) >> 12;  //rest delay
				if (a) {
					block[z][y][x] = block[z][y][x] & 0x03cfff;  //clear rest delay 
					block[z][y][x] += ((a - 1) << 12) + ((upbs == 3) ? 0x008000 : 0x020000);
					addtolist(x, y, z);
				}
				else {
					b = block[z][y][x] & 0x010000;  //level next
					block[z][y][x] = (block[z][y][x] & 0x02bfff) + (b >> 2);
				}
			}
		}
		if ((block[z][y][x] & ((upbs == 3) ? 0x020000 : 0x008000)) == 0 && ((block[z][y][x] & ((upbs == 3) ? 0x008000 : 0x020000)) == 0 || (block[z][y][x] & 0x3000) >> 2 == (block[z][y][x] & 0x0c00))) {
			f = (block[z][y][x] & 0x0300) >> 8;
			switch (f) {
			case 0:
				x1 = x - 1;
				z1 = z;
				break;
			case 1:
				x1 = x + 1;
				z1 = z;
				break;
			case 2:
				x1 = x;
				z1 = z - 1;
				break;
			case 3:
				x1 = x;
				z1 = z + 1;
				break;
			}
			a = glevdio(x1, y, z1, f, block);
			if (a && (block[z][y][x] & 0x4000) == 0) {
				block[z][y][x] = (block[z][y][x] & 0x4fff) + ((upbs == 3) ? 0x008000 : 0x020000) + 0x010000 + ((block[z][y][x] & 0x0c00) << 2);
				addtolist(x, y, z);
			}
			else if (a == 0 && (block[z][y][x] & 0x4000)) {
				block[z][y][x] = (block[z][y][x] & 0x4fff) + ((upbs == 3) ? 0x008000 : 0x020000) + ((block[z][y][x] & 0x0c00) << 2);
				addtolist(x, y, z);
			}
			else {
				block[z][y][x] = (block[z][y][x] & 0x4fff);
				removefromlist(x, y, z);
			}
		}
	}
	if ((ib & 0x4fff) != (block[z][y][x] & 0x4fff) || s) {
		regionsupdate[x / 16 + z / 16 * 16].opaque |= 1;
		if (x > 0) {
			upb(x - 1, y, z, 0, block, regionsupdate);
		}
		if (y > 0) {
			upb(x, y - 1, z, 0, block, regionsupdate);
		}
		if (z > 0) {
			upb(x, y, z - 1, 0, block, regionsupdate);
		}
		if (x < 255) {
			upb(x + 1, y, z, 0, block, regionsupdate);
		}
		if (y < 255) {
			upb(x, y + 1, z, 0, block, regionsupdate);
		}
		if (z < 255) {
			upb(x, y, z + 1, 0, block, regionsupdate);
		}
		return 1;
	}
	return 0;
}

int uprt(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate) {
	int ib, x1, y1, z1, i;
	ib = block[z][y][x];
	if (block[z][y][x] & ((upbs == 3) ? 0x4000 : 0x1000)) {
		block[z][y][x] = (block[z][y][x] & 0xffff87ff) + ((ib & 0x2000) >> 2);
	}
	x1 = x;
	z1 = z;
	y1 = y;
	switch ((block[z][y][x] & 0x0700) >> 8) {
	case 0:
		y1 -= 1;
		break;
	case  1:
		x1 += 1;
		break;
	case 2:
		x1 -= 1;
		break;
	case 3:
		z1 += 1;
		break;
	case 4:
		z1 -= 1;
		break;
	}
	if (y1 < 0 || tr[block[z1][y1][x1] & 0x00ff]) {
		block[z][y][x] = 0;
	}
	else {
		if ((block[z][y][x] & 0x0800) && (block[z1][y1][x1] & 0x0300)) {
			block[z][y][x] = (block[z][y][x] & 0x0fff) + ((upbs == 3) ? 0x1000 : 0x4000);
			addtolist(x, y, z);
		}
		else if ((block[z][y][x] & 0x0800) == 0 && (block[z1][y1][x1] & 0x0300) == 0) {
			block[z][y][x] = (block[z][y][x] & 0x0fff) + ((upbs == 3) ? 0x3000 : 0x6000);
			addtolist(x, y, z);
		}
		else {
			block[z][y][x] &= 0x0fff;
			removefromlist(x, y, z);
		}
	}
	if ((ib & 0x0fff) != (block[z][y][x] & 0x0fff) || s) {
		regionsupdate[x / 16 + z / 16 * 16].opaque |= 1;
		if (x > 0) {
			upb(x - 1, y, z, 0, block, regionsupdate);
		}
		if (y > 0) {
			upb(x, y - 1, z, 0, block, regionsupdate);
		}
		if (z > 0) {
			upb(x, y, z - 1, 0, block, regionsupdate);
		}
		if (x < 255) {
			upb(x + 1, y, z, 0, block, regionsupdate);
		}
		if (y < 255) {
			upb(x, y + 1, z, 0, block, regionsupdate);
		}
		if (z < 255) {
			upb(x, y, z + 1, 0, block, regionsupdate);
		}
		return 1;
	}
	return 0;
}

int upitb(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate) {
	int ib;
	ib = block[z][y][x];
	block[z][y][x] &= 0xfffffcff;
	if ((block[z][y][x] & 0x00ff) == 3) {
		block[z][y][x] |= 0x0200;
	}
	if (y < 255 && (block[z][y + 1][x] & 0x00ff) == 4 && (block[z][y + 1][x] & 0xf000)) {
		block[z][y][x] |= 0x0100;
	}
	if (y>0 && (block[z][y - 1][x] & 0x00ff) == 5 && (block[z][y - 1][x] & 0x0800)) {
		block[z][y][x] |= 0x0200;
	}
	if (x > 0 && (block[z][y][x - 1] & 0x0fff) == 0x0404 && (block[z][y][x - 1] & 0xf000)) {
		block[z][y][x] |= 0x0100;
	}
	if (z > 0 && (block[z - 1][y][x] & 0x0fff) == 0x0804 && (block[z - 1][y][x] & 0xf000)) {
		block[z][y][x] |= 0x0100;
	}
	if (x < 255 && (block[z][y][x + 1] & 0x0fff) == 0x0104 && (block[z][y][x + 1] & 0xf000)) {
		block[z][y][x] |= 0x0100;
	}
	if (z < 255 && (block[z + 1][y][x] & 0x0fff) == 0x0204 && (block[z + 1][y][x] & 0xf000)) {
		block[z][y][x] |= 0x0100;
	}
	if (x > 0 && (block[z][y][x - 1] & 0x03ff) == 0x0006 && (block[z][y][x - 1] & 0x4000)) {
		block[z][y][x] |= 0x0200;
	}
	if (z > 0 && (block[z - 1][y][x] & 0x03ff) == 0x0206 && (block[z - 1][y][x] & 0x4000)) {
		block[z][y][x] |= 0x0200;
	}
	if (x < 255 && (block[z][y][x + 1] & 0x03ff) == 0x0106 && (block[z][y][x + 1] & 0x4000)) {
		block[z][y][x] |= 0x0200;
	}
	if (z < 255 && (block[z + 1][y][x] & 0x03ff) == 0x0306 && (block[z + 1][y][x] & 0x4000)) {
		block[z][y][x] |= 0x0200;
	}
	if (s || block[z][y][x] != ib) {
		if (x > 0)
			upb(x - 1, y, z, 0, block, regionsupdate);
		if (y > 0)
			upb(x, y - 1, z, 0, block, regionsupdate);
		if (z > 0)
			upb(x, y, z - 1, 0, block, regionsupdate);
		if (x < 255)
			upb(x + 1, y, z, 0, block, regionsupdate);
		if (y < 255)
			upb(x, y + 1, z, 0, block, regionsupdate);
		if (z < 255)
			upb(x, y, z + 1, 0, block, regionsupdate);
	}
	return 0;
}

int uprl(int x, int y, int z, int s, int(*block)[256][256], elementtypes *regionsupdate) {
	int ir, lev;
	ir = block[z][y][x];
	lev = 0;
	if (y == 0) {
		block[z][y][x] = 0;
	}
	else if (tr[block[z][y - 1][x] & 0x00ff]) {
		block[z][y][x] = 0;
	}
	else {
		block[z][y][x] &= 0xffff00ff;
		if (x < 255) {
			lev = glev(x + 1, y, z, 1, block);
			if ((block[z][y][x + 1] & 0x00ff) == 4 || (block[z][y][x + 1] & 0x00ff) == 5) {
				block[z][y][x] |= 0x0100;
			}
			else if (tr[block[z][y][x + 1] & 0x00ff] && y > 0 && (block[z][y - 1][x + 1] & 0x00ff) == 4) {
				block[z][y][x] |= 0x0100;
				lev = glev(x + 1, y - 1, z, 1, block);
			}
			else if (!tr[block[z][y][x + 1] & 0x00ff] && y < 255 && (block[z][y + 1][x + 1] & 0x00ff) == 4 && tr[block[z][y + 1][x] & 0x00ff]) {
				block[z][y][x] |= 0x0100;
				lev = glev(x + 1, y + 1, z, 1, block);
			}
		}
		if (z < 255) {
			lev = max(lev, glev(x, y, z + 1, 3, block));
			if ((block[z + 1][y][x] & 0x00ff) == 4 || (block[z + 1][y][x] & 0x00ff) == 5) {
				block[z][y][x] |= 0x0200;
			}
			else if (tr[block[z + 1][y][x] & 0x00ff] && y > 0 && (block[z + 1][y - 1][x] & 0x00ff) == 4) {
				block[z][y][x] |= 0x0200;
				lev = max(lev, glev(x, y - 1, z + 1, 3, block));
			}
			else if (!tr[block[z + 1][y][x] & 0x00ff] && y < 255 && (block[z + 1][y + 1][x] & 0x00ff) == 4 && tr[block[z][y + 1][x] & 0x00ff]) {
				block[z][y][x] |= 0x0200;
				lev = max(lev, glev(x, y + 1, z + 1, 3, block));
			}
		}
		if (x > 0) {
			lev = max(lev, glev(x - 1, y, z, 0, block));
			if ((block[z][y][x - 1] & 0x00ff) == 4 || (block[z][y][x - 1] & 0x00ff) == 5) {
				block[z][y][x] |= 0x0400;
			}
			else if (tr[block[z][y][x - 1] & 0x00ff] && y > 0 && (block[z][y - 1][x - 1] & 0x00ff) == 4) {
				block[z][y][x] |= 0x0400;
				lev = max(lev, glev(x - 1, y - 1, z, 0, block));
			}
			else if (!tr[block[z][y][x - 1] & 0x00ff] && y < 255 && (block[z][y + 1][x - 1] & 0x00ff) == 4 && tr[block[z][y + 1][x] & 0x00ff]) {
				block[z][y][x] |= 0x0400;
				lev = max(lev, glev(x - 1, y + 1, z, 0, block));
			}
		}
		if (z > 0) {
			lev = max(lev, glev(x, y, z - 1, 2, block));
			if ((block[z - 1][y][x] & 0x00ff) == 4 || (block[z - 1][y][x] & 0x00ff) == 5) {
				block[z][y][x] |= 0x0800;
			}
			else if (tr[block[z - 1][y][x] & 0x00ff] && y > 0 && (block[z - 1][y - 1][x] & 0x00ff) == 4) {
				block[z][y][x] |= 0x0800;
				lev = max(lev, glev(x, y - 1, z - 1, 2, block));
			}
			else if (!tr[block[z - 1][y][x] & 0x00ff] && y < 255 && (block[z - 1][y + 1][x] & 0x00ff) == 4 && tr[block[z][y + 1][x] & 0x00ff]) {
				block[z][y][x] |= 0x0800;
				lev = max(lev, glev(x, y + 1, z - 1, 2, block));
			}
		}
		lev = max(lev, max(glev(x, y + 1, z, 4, block), glev(x, y - 1, z, 5, block)));
		block[z][y][x] += lev << 12;
	}
	if (ir != block[z][y][x] || s) {
		regionsupdate[x / 16 + z / 16 * 16].opaque |= 1;
		if (x > 0) {
			upb(x - 1, y, z, 0, block, regionsupdate);
			if (y > 0) {
				upb(x - 1, y - 1, z, 0, block, regionsupdate);
			}
			if (y < 255) {
				upb(x - 1, y + 1, z, 0, block, regionsupdate);
			}
		}
		if (y > 0) {
			upb(x, y - 1, z, 0, block, regionsupdate);
		}
		if (z > 0) {
			upb(x, y, z - 1, 0, block, regionsupdate);
			if (y > 0) {
				upb(x, y - 1, z - 1, 0, block, regionsupdate);
			}
			if (y < 255) {
				upb(x, y + 1, z - 1, 0, block, regionsupdate);
			}
		}
		if (x < 255) {
			upb(x + 1, y, z, 0, block, regionsupdate);
			if (y > 0) {
				upb(x + 1, y - 1, z, 0, block, regionsupdate);
			}
			if (y < 255) {
				upb(x + 1, y + 1, z, 0, block, regionsupdate);
			}
		}
		if (y < 255) {
			upb(x, y + 1, z, 0, block, regionsupdate);
		}
		if (z < 255) {
			upb(x, y, z + 1, 0, block, regionsupdate);
			if (y > 0) {
				upb(x, y - 1, z + 1, 0, block, regionsupdate);
			}
			if (y < 255) {
				upb(x, y + 1, z + 1, 0, block, regionsupdate);
			}
		}
		return 1;
	}
	return 0;
}

int glev(int x, int y, int z, int f, int(*block)[256][256]) {
	int n;
	if (x < 0 || x>255 || z < 0 || z>255 || y < 0 || y>255) {
		return 0;
	}
	n = block[z][y][x] & 0x00ff;
	switch (n) {
	case 0:
		return 0;
	case 4:
		return max(0, ((block[z][y][x] & 0xf000) >> 12) - 1);
	case 5:
		return (block[z][y][x] & 0x0800) ? 15 : 0;
	case 6:
		return (((block[z][y][x] & 0x0300) >> 8) == f && (block[z][y][x] & 0x4000)) ? 15 : 0;
	default:
		return (block[z][y][x] & 0x0200) ? 15 : 0;
	}
}

int glevdio(int x, int y, int z, int f, int(*block)[256][256]) {
	int n;
	if (x < 0 || x>255 || z < 0 || z>255 || y < 0 || y>255) {
		return 0;
	}
	n = block[z][y][x] & 0x00ff;
	switch (n) {
	case 0:
		return 0;
	case 4:
		return (block[z][y][x] & 0xf000) >> 12;
	case 5:
		return (block[z][y][x] & 0x0800) ? 15 : 0;
	case 6:
		return (((block[z][y][x] & 0x0300) >> 8) == f && (block[z][y][x] & 0x4000)) ? 15 : 0;
	default:
		return (block[z][y][x] & 0x0300) ? 15 : 0;
	}
}

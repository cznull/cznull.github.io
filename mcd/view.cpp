#include "stdafx.h"

struct buffer3 {
	int max;
	int count;
	float3 *vertex;
	float2 *texture;
	int vertexcount[256];
};

buffer3 vb;

int whichsurface(movestate *pe, int(*block)[256][256]);

inline void rotatef(float3 axis, float3 &pos, float angle) {
	float a, cosang, sinang;
	float3 rad1, rad2, axi;
	a = axis.x * pos.x + axis.y*pos.y + axis.z * pos.z;
	axi = { a*axis.x,a*axis.y,a*axis.z };
	rad1 = { pos.x - axi.x,pos.y - axi.y ,pos.z - axi.z };
	rad2 = { axis.y*rad1.z - axis.z*rad1.y, axis.z*rad1.x - axis.x*rad1.z, axis.x*rad1.y - axis.y*rad1.x };
	cosang = cosf(angle);
	sinang = sinf(angle);
	pos = { axi.x + cosang * rad1.x + sinang * rad2.x,axi.y + cosang * rad1.y + sinang * rad2.y ,axi.z + cosang * rad1.z + sinang * rad2.z };
}

const int texsize = 64;
void(*show[256])(int, int, int, int);

const float texms[100] = { 0.0,0.25,0.25,0,0.25,0.50,0.50,0.25,0,0.25,0.25,0 , 0,0.25,0.25,0,0,0.25,0.25,0,0,0.25,0.25,0 , 0.25,0.50,0.50,0.25,0.25,0.50,0.50,0.25,0.25,0.50,0.50,0.25 };
const float texmt[100] = { 0.50,0.50,0.25,0.25, 0.25,0.25,0,0,0.25,0.25,0,0 , 0.25,0.25,0,0, 0.25,0.25,0,0,0.25,0.25,0,0 , 0.50,0.50,0.25,0.25,0.50,0.50,0.25,0.25,0.50,0.50,0.25,0.25 };
const float rverx[64] = { 0.375,0.625,0.625,0.375, 0.000,1.000,1.000,0.000, 0.375,0.625,0.625,0.375, 0.375,1.000,1.000,0.375, 0.000,1.000,1.000,0.000, 0.000,1.000,1.000,0.000, 0.000,0.625,0.625,0.000, 0.000,1.000,1.000,0.000, 0.375,0.625,0.625,0.375, 0.375,1.000,1.000,0.375, 0.375,0.625,0.625,0.375, 0.375,1.000,1.000,0.375, 0.000,0.625,0.625,0.000, 0.000,1.000,1.000,0.000, 0.000,0.625,0.625,0.000, 0.000,1.000,1.000,0.000 };
const float rverz[64] = { 0.375,0.375,0.625,0.625, 0.375,0.375,0.625,0.625, 0.000,0.000,1.000,1.000, 0.375,0.375,1.000,1.000, 0.375,0.375,0.625,0.625, 0.375,0.375,0.625,0.625, 0.375,0.375,1.000,1.000, 0.375,0.375,1.000,1.000, 0.000,0.000,1.000,1.000, 0.000,0.000,0.625,0.625, 0.000,0.000,1.000,1.000, 0.000,0.000,1.000,1.000, 0.000,0.000,0.625,0.625, 0.000,0.000,0.625,0.625, 0.000,0.000,1.000,1.000, 0.000,0.000,1.000,1.000 };
const float rtexms[8] = { 0.59375,0.65625,0.65625,0.59375,0.609375,0.640625,0.640625,0.609375 };
const float rtexmt[8] = { 0.25,0.25,0.078125,0.078125,0.125,0.125,0.09375,0.09375 };

int whichsurface(movestate *pe, int(*block)[256][256]) {
	if (tr[block[pe->point.blockz][pe->point.blocky][pe->point.blockx] & 0x00ff] == 0) {
		pe->point.ispoint = 1;
		return 1;
	}
	switch (block[pe->point.blockz][pe->point.blocky][pe->point.blockx] & 0x00ff) {
	case 4:
	case 6:
		if (pe->point.sury == pe->point.blocky - 1) {
			pe->point.ispoint = 2;
			return 1;
		}
		if (pe->point.surx == pe->point.blockx - 1 && pe->at.y*(pe->point.blockx - pe->pos.x) / pe->at.x + pe->pos.y <= pe->point.blocky + 0.125) {
			pe->point.ispoint = 2;
			return 1;
		}
		if (pe->point.surx == pe->point.blockx + 1 && pe->at.y*(pe->point.blockx + 1 - pe->pos.x) / pe->at.x + pe->pos.y <= pe->point.blocky + 0.125) {
			pe->point.ispoint = 2;
			return 1;
		}
		if (pe->point.surz == pe->point.blockz - 1 && pe->at.y*(pe->point.blockz - pe->pos.z) / pe->at.z + pe->pos.y <= pe->point.blocky + 0.125) {
			pe->point.ispoint = 2;
			return 1;
		}
		if (pe->point.surz == pe->point.blockz + 1 && pe->at.y*(pe->point.blockz + 1 - pe->pos.z) / pe->at.z + pe->pos.y <= pe->point.blocky + 0.125) {
			pe->point.ispoint = 2;
			return 1;
		}
		if (pe->at.y < 0 && pe->at.x*(pe->point.blocky + 0.125 - pe->pos.y) / pe->at.y + pe->pos.x <= pe->point.blockx + 1 && pe->at.x*(pe->point.blocky + 0.125 - pe->pos.y) / pe->at.y + pe->pos.x >= pe->point.blockx && pe->at.z*(pe->point.blocky + 0.125 - pe->pos.y) / pe->at.y + pe->pos.z <= pe->point.blockz + 1 && pe->at.z*(pe->point.blocky + 0.125 - pe->pos.y) / pe->at.y + pe->pos.z >= pe->point.blockz) {
			pe->point.sury = pe->point.blocky + 1;
			pe->point.surx = pe->point.blockx;
			pe->point.surz = pe->point.blockz;
			pe->point.ispoint = 2;
			return 1;
		}
		return 0;
	case 5:
		pe->point.ispoint = 1;
		return 1;
	}
}

void updatepoint(movestate *pe, int(*block)[256][256]) {
	int xor0, yor0, zor0;
	double dx, dy, dz;
	pe->point.blockx = (int)pe->pos.x;
	pe->point.blocky = (int)pe->pos.y;
	pe->point.blockz = (int)pe->pos.z;
	xor0 = pe->at.x > 0 ? 1 : 0;
	yor0 = pe->at.y > 0 ? 1 : 0;
	zor0 = pe->at.z > 0 ? 1 : 0;
	pe->point.ispoint = 0;
	while ((((xor0) && (pe->point.blockx <= 255)) || ((!xor0) && (pe->point.blockx >= 0))) && (((yor0) && (pe->point.blocky <= 255)) || ((!yor0) && (pe->point.blocky >= 0))) && (((zor0) && (pe->point.blockz <= 255)) || ((!zor0) && (pe->point.blockz >= 0)))) {
		pe->point.surx = pe->point.blockx;
		pe->point.sury = pe->point.blocky;
		pe->point.surz = pe->point.blockz;
		dx = double(pe->point.blockx + xor0) - pe->pos.x;
		dy = double(pe->point.blocky + yor0) - pe->pos.y;
		dz = double(pe->point.blockz + zor0) - pe->pos.z;
		if (dx*pe->at.y*dx*pe->at.y < dy*pe->at.x*dy*pe->at.x) {
			if (dx*dx*pe->at.z*pe->at.z < dz*dz*pe->at.x*pe->at.x)
				pe->point.blockx += (pe->at.x > 0) ? 1 : -1;
			else
				pe->point.blockz += (pe->at.z > 0) ? 1 : -1;
		}
		else {
			if (dy*dy*pe->at.z*pe->at.z < dz*dz*pe->at.y*pe->at.y)
				pe->point.blocky += (pe->at.y > 0) ? 1 : -1;
			else
				pe->point.blockz += (pe->at.z > 0) ? 1 : -1;
		}
		if (0 <= pe->point.blockx && pe->point.blockx <= 255 && 0 <= pe->point.blocky && pe->point.blocky <= 255 && 0 <= pe->point.blockz && pe->point.blockz <= 255 && block[pe->point.blockz][pe->point.blocky][pe->point.blockx]) {
			if (whichsurface(pe, block)) {
				break;
			}
		}
	}
}

void addpoint(float x, float y, float z, float s, float t) {
	if (vb.count >= vb.max) {
		float3 *v1;
		float2 *t1;
		v1 = (float3*)malloc(vb.max * sizeof(float3) * 2);
		memcpy(v1, vb.vertex, vb.max * sizeof(float3));
		free(vb.vertex);
		vb.vertex = v1;
		t1 = (float2*)malloc(vb.max * sizeof(float2) * 2);
		memcpy(t1, vb.texture, vb.max * sizeof(float2));
		free(vb.texture);
		vb.texture = t1;
		vb.max *= 2;
	}
	vb.vertex[vb.count].x = x;
	vb.vertex[vb.count].y = y;
	vb.vertex[vb.count].z = z;
	vb.texture[vb.count].s = s;
	vb.texture[vb.count].t = t;
	vb.count++;
}

void showblock(int k, int j, int i, int n, int(*block)[256][256], int buffernumber) {
	int blockid;
	blockid = block[i][j][k];
	switch (blockid & 0xff) {
	case 1:
	case 2:
	case 3:
		if (k < 255 && tr[block[i][j][k + 1] & 0x00ff] || k == 255) {
			addpoint(k + 1, j, i + 1, texms[n * 12 - 12 + 0], texmt[n * 12 - 12 + 0]);
			addpoint(k + 1, j, i, texms[n * 12 - 12 + 1], texmt[n * 12 - 12 + 1]);
			addpoint(k + 1, j + 1, i, texms[n * 12 - 12 + 2], texmt[n * 12 - 12 + 2]);
			addpoint(k + 1, j + 1, i + 1, texms[n * 12 - 12 + 3], texmt[n * 12 - 12 + 3]);
		}
		if (k > 0 && tr[block[i][j][k - 1] & 0x00ff] || k == 0) {
			addpoint(k, j, i, texms[n * 12 - 12 + 0], texmt[n * 12 - 12 + 0]);
			addpoint(k, j, i + 1, texms[n * 12 - 12 + 1], texmt[n * 12 - 12 + 1]);
			addpoint(k, j + 1, i + 1, texms[n * 12 - 12 + 2], texmt[n * 12 - 12 + 2]);
			addpoint(k, j + 1, i, texms[n * 12 - 12 + 3], texmt[n * 12 - 12 + 3]);
		}
		if (j < 255 && tr[block[i][j + 1][k] & 0x00ff] || j == 255) {
			addpoint(k + 1, j + 1, i + 1, texms[n * 12 - 12 + 4], texmt[n * 12 - 12 + 4]);
			addpoint(k + 1, j + 1, i, texms[n * 12 - 12 + 5], texmt[n * 12 - 12 + 5]);
			addpoint(k, j + 1, i, texms[n * 12 - 12 + 6], texmt[n * 12 - 12 + 6]);
			addpoint(k, j + 1, i + 1, texms[n * 12 - 12 + 7], texmt[n * 12 - 12 + 7]);
		}
		if (j > 0 && tr[block[i][j - 1][k] & 0x00ff] || j == 0) {
			addpoint(k, j, i, texms[n * 12 - 12 + 8], texmt[n * 12 - 12 + 8]);
			addpoint(k + 1, j, i, texms[n * 12 - 12 + 9], texmt[n * 12 - 12 + 9]);
			addpoint(k + 1, j, i + 1, texms[n * 12 - 12 + 10], texmt[n * 12 - 12 + 10]);
			addpoint(k, j, i + 1, texms[n * 12 - 12 + 11], texmt[n * 12 - 12 + 11]);
		}
		if (i < 255 && tr[block[i + 1][j][k] & 0x00ff] || i == 255) {
			addpoint(k, j, i + 1, texms[n * 12 - 12 + 0], texmt[n * 12 - 12 + 0]);
			addpoint(k + 1, j, i + 1, texms[n * 12 - 12 + 1], texmt[n * 12 - 12 + 1]);
			addpoint(k + 1, j + 1, i + 1, texms[n * 12 - 12 + 2], texmt[n * 12 - 12 + 2]);
			addpoint(k, j + 1, i + 1, texms[n * 12 - 12 + 3], texmt[n * 12 - 12 + 3]);
		}
		if (i > 0 && tr[block[i - 1][j][k] & 0x00ff] || i == 0) {
			addpoint(k + 1, j, i, texms[n * 12 - 12 + 0], texmt[n * 12 - 12 + 0]);
			addpoint(k, j, i, texms[n * 12 - 12 + 1], texmt[n * 12 - 12 + 1]);
			addpoint(k, j + 1, i, texms[n * 12 - 12 + 2], texmt[n * 12 - 12 + 2]);
			addpoint(k + 1, j + 1, i, texms[n * 12 - 12 + 3], texmt[n * 12 - 12 + 3]);
		}
		break;
	case 4: {
		int bright, connect;
		bright = (blockid & 0xF000) >> 12;
		connect = (blockid & 0x0F00) >> 8;
		addpoint(k + rverx[connect * 4 + 0], j + 0.03, i + rverz[connect * 4 + 0], (bright + 0.5) / 64, 32.5 / 64);
		addpoint(k + rverx[connect * 4 + 1], j + 0.03, i + rverz[connect * 4 + 1], (bright + 0.5) / 64, 32.5 / 64);
		addpoint(k + rverx[connect * 4 + 2], j + 0.03, i + rverz[connect * 4 + 2], (bright + 0.5) / 64, 32.5 / 64);
		addpoint(k + rverx[connect * 4 + 3], j + 0.03, i + rverz[connect * 4 + 3], (bright + 0.5) / 64, 32.5 / 64);
		break;
	}
	case 5: {
		int a, b, c, m;
		float3 ver[20] = { { -0.125,-1.0,0.0625 },{ 0.125,-1.0,0.0625 },{ 0.125,-0.3125,0.0625 },{ -0.125,-0.3125,0.0625 },{ 0.125,-1.0,-0.0625 },{ -0.125,-1.0,-0.0625 },{ -0.125,-0.3125,-0.0625 },{ 0.125,-0.3125,-0.0625 },
						   { 0.0625,-1.0,0.125 },{ 0.0625,-1.0,-0.125 },{ 0.0625,-0.3125,-0.125 },{ 0.0625,-0.3125,0.125 },{ -0.0625,-1.0,-0.125 },{ -0.0625,-1.0,0.125 },{ -0.0625,-0.3125,0.125 },{ -0.0625,-0.3125,-0.125 },
						   { 0.0625,-0.375,0.0625 },{ 0.0625,-0.375,-0.0625 },{ -0.0625,-0.375,-0.0625 },{ -0.0625,-0.375,0.0625 } };
		a = (block[i][j][k] & 0x0700) >> 8;
		b = (block[i][j][k] & 0x0800) >> 11;
		switch (a) {
		case 0:
			break;
		case 1:
			for (m = 0; m < 20; m++) {
				rotatef({ 0.0,0.0,1.0 }, ver[m], PI / 6);
			}
			break;
		case 2:
			for (m = 0; m < 20; m++) {
				rotatef({ 0.0,0.0,-1.0 }, ver[m], PI / 6);
			}
			break;
		case 3:
			for (m = 0; m < 20; m++) {
				rotatef({ -1.0,0.0,0.0 }, ver[m], PI / 6);
			}
			break;
		case 4:
			for (m = 0; m < 20; m++) {
				rotatef({ 1.0,0.0,0.0 }, ver[m], PI / 6);
			}
			break;
		}
		for (m = 0; m < 16; m++) {
			addpoint(ver[m].x + k + 0.5, ver[m].y + j + 1.0, ver[m].z + i + 0.5, rtexms[m % 4], rtexmt[m % 4] + b * 0.25);
		}
		for (; m < 20; m++) {
			addpoint(ver[m].x + k + 0.5, ver[m].y + j + 1.0, ver[m].z + i + 0.5, rtexms[m % 4 + 4], rtexmt[m % 4 + 4] + b * 0.25);
		}
		break;
	}
	case 6: {
		int a, b;
		float3 ver[4] = { {-0.5,0,-0.5},{-0.5,0,0.5}, {0.5,0,0.5}, {0.5,0,-0.5} };
		a = (block[i][j][k] & 0x0300) >> 8;
		b = (block[i][j][k] & 0x4000) >> 14;
		switch (a) {
		case 0:
			break;
		case 1:
			rotatef({ 0.0,1.0,0.0 }, ver[0], PI);
			rotatef({ 0.0,1.0,0.0 }, ver[1], PI);
			rotatef({ 0.0,1.0,0.0 }, ver[2], PI);
			rotatef({ 0.0,1.0,0.0 }, ver[3], PI);
			break;
		case 2:
			rotatef({ 0.0,1.0,0.0 }, ver[0], PI*1.5);
			rotatef({ 0.0,1.0,0.0 }, ver[1], PI*1.5);
			rotatef({ 0.0,1.0,0.0 }, ver[2], PI*1.5);
			rotatef({ 0.0,1.0,0.0 }, ver[3], PI*1.5);
			break;
		case 3:
			rotatef({ 0.0,1.0,0.0 }, ver[0], PI*0.5);
			rotatef({ 0.0,1.0,0.0 }, ver[1], PI*0.5);
			rotatef({ 0.0,1.0,0.0 }, ver[2], PI*0.5);
			rotatef({ 0.0,1.0,0.0 }, ver[3], PI*0.5);
			break;
		}
		addpoint(ver[0].x + k + 0.5, ver[0].y + j + 0.03125, ver[0].z + i + 0.5, 0.75, 0.25);
		addpoint(ver[1].x + k + 0.5, ver[1].y + j + 0.03125, ver[1].z + i + 0.5, 1.00, 0.25);
		addpoint(ver[2].x + k + 0.5, ver[2].y + j + 0.03125, ver[2].z + i + 0.5, 1.00, 0.00);
		addpoint(ver[3].x + k + 0.5, ver[3].y + j + 0.03125, ver[3].z + i + 0.5, 0.75, 0.00);
		break;
	}
	}
}

void upregionsurface(int(*block)[256][256], int rx, int rz, int t, GLuint *vbo, GLuint *tbo) {
	int i, j, k, n, buffernumber;
	if (rx < 0 || rz < 0 || rx>16 || rz>15) {
		return;
	}
	buffernumber = rz * 16 + rx;
	vb.count = 0;
	for (i = rz * 16; i < rz * 16 + 16; i++) {
		for (j = 0; j < 256; j++) {
			for (k = rx * 16; k < rx * 16 + 16; k++) {
				n = (block[i][j][k]) % 256;
				//(*show[n])(k, j, i, block[i][j][k] & 0x00ff);
				if (n) {
					showblock(k, j, i, block[i][j][k] & 0x00ff, block, buffernumber);
				}
			}
		}
	}
	vb.vertexcount[buffernumber] = vb.count;
	glBindBuffer(GL_ARRAY_BUFFER, vbo[buffernumber]);
	glBufferData(GL_ARRAY_BUFFER, vb.count * sizeof(float3), vb.vertex, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, tbo[buffernumber]);
	glBufferData(GL_ARRAY_BUFFER, vb.count * sizeof(float2), vb.texture, GL_STATIC_DRAW);
}

int initview(GLuint *tex, GLuint *vbo, GLuint *tbo, int(*block)[256][256], HBITMAP teximg) {
	int i, j;
	ubyte4 *texdata;
	glGenBuffers(256, vbo);
	glGenBuffers(256, tbo);
	glGenTextures(1, tex);

	glBindTexture(GL_TEXTURE_2D, *tex);
	texdata = (ubyte4*)malloc(texsize * texsize * sizeof(ubyte4));
	GetBitmapBits(teximg, texsize*texsize*sizeof(ubyte4), texdata);
	for (i = 0; i < texsize*texsize; i++) {
		if (texdata[i].x != 255 || texdata[i].z != 255 || texdata[i].z != 255) {
			texdata[i].w = 255;
		}
		else {
			texdata[i].w = 0;
		}
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texsize, texsize, 0, GL_BGRA, GL_UNSIGNED_BYTE, texdata);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	free(texdata);

	vb.vertex = (float3*)malloc(128 * sizeof(float3));
	vb.texture = (float2*)malloc(128 * sizeof(float2));
	vb.max = 128;

	for (i = 0; i < 256; i++) {
		upregionsurface(block, i % 16, i / 16, 0, vbo, tbo);
	}
	return 1;
}

int draw(movestate *pe,GLuint *vbo, GLuint *tbo) {
	int i;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(pe->pos.x, pe->pos.y - 0.1*(pe->onlandstate.isonland && (pe->up == -1)), pe->pos.z, cos(pe->ang2)*cos(pe->ang1) + pe->pos.x, sin(pe->ang2) + pe->pos.y - 0.1*(pe->onlandstate.isonland && (pe->up == -1)), sin(pe->ang1)*cos(pe->ang2) + pe->pos.z, 0, cos(pe->ang2), 0);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(0x00004100);
	glColor3f(1.0f, 1.0f, 1.0f);
	for (i = 0; i < 256; i++) {
		glBindBuffer(GL_ARRAY_BUFFER, vbo[i]);
		glVertexPointer(3, GL_FLOAT, 0, NULL);
		glBindBuffer(GL_ARRAY_BUFFER, tbo[i]);
		glTexCoordPointer(2, GL_FLOAT, 0, NULL);
		glDrawArrays(GL_QUADS, 0, vb.vertexcount[i]);
	}
	glLoadIdentity();
	glScalef(0.995, 0.995, 0.995);
	gluLookAt(pe->pos.x, pe->pos.y - 0.1*(pe->onlandstate.isonland && (pe->up == -1)), pe->pos.z, cos(pe->ang2)*cos(pe->ang1) + pe->pos.x, sin(pe->ang2) + pe->pos.y - 0.1*(pe->onlandstate.isonland && (pe->up == -1)), sin(pe->ang1)*cos(pe->ang2) + pe->pos.z, 0, cos(pe->ang2), 0);
	switch (pe->point.ispoint) {
	case 1:
		glTranslatef(pe->point.blockx, pe->point.blocky, pe->point.blockz);
		glBegin(GL_LINE_LOOP);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(-0.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, 1.00f, 1.00f);
		glVertex3f(1.00f, 1.00f, 1.00f);
		glVertex3f(1.00f, 1.00f, -0.00f);
		glVertex3f(-0.00f, 1.00f, -0.00f);
		glEnd();
		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, 1.00f, -0.00f);
		glVertex3f(1.00f, -0.00f, 1.00f);
		glVertex3f(1.00f, 1.00f, 1.00f);
		glVertex3f(-0.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, -0.00f, -0.00f);
		glVertex3f(-0.00f, 1.00f, 1.00f);
		glVertex3f(-0.00f, 1.00f, -0.00f);
		glEnd();
		break;
	case 2:
		glTranslatef(pe->point.blockx, pe->point.blocky, pe->point.blockz);
		glBegin(GL_LINE_LOOP);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(-0.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, 0.125f, 1.00f);
		glVertex3f(1.00f, 0.125f, 1.00f);
		glVertex3f(1.00f, 0.125f, -0.00f);
		glVertex3f(-0.00f, 0.125f, -0.00f);
		glEnd();
		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.00f, -0.00f, -0.00f);
		glVertex3f(1.00f, 0.125f, -0.00f);
		glVertex3f(1.00f, -0.00f, 1.00f);
		glVertex3f(1.00f, 0.125f, 1.00f);
		glVertex3f(-0.00f, -0.00f, 1.00f);
		glVertex3f(-0.00f, -0.00f, -0.00f);
		glVertex3f(-0.00f, 0.125f, 1.00f);
		glVertex3f(-0.00f, 0.125f, -0.00f);
		glEnd();
		break;
	}
	glLoadIdentity();
	glClear(0x00000100);
	glBegin(GL_LINES);
	glColor3f(0, 0, 0);
	glVertex3f(0.2f, 0.0f, -10.0f);
	glVertex3f(-0.2f, 0.0f, -10.0f);
	glVertex3f(0.0f, 0.2f, -10.0f);
	glVertex3f(0.0f, -0.2f, -10.0f);
	glEnd();
	SwapBuffers(wglGetCurrentDC());
	return 1;
}
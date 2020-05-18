#include <iostream>
#include <stdio.h>
#include <io.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <spine/spine.h>
#include <lz4.h>
#pragma comment(lib,"liblz4_static.lib")

struct block_t {
	int headpos;
	int uncom, com, flag;
};

struct file_t {
	long long length, offset;
	std::string name;
	int flag;
};

struct typeTreeNode_t {
	int Version, Level, IsArray, TypeStrOffset, NameStrOffset, ByteSize, Index, MetaFlag;
	std::string Type, Name;
};


struct type_t {
	char ScriptID[16], OldTypeHash[16];
	int classID, IsStrippedType, ScriptTypeIndex;
	std::vector<typeTreeNode_t> typeTree;
};

struct objectInfo_t {
	long long PathID;
	int byteStart, byteSize, typeID;
	type_t type;
};

struct imgname_t {
	std::string rgb, a, name;
	long long rgbid, aid;
};

struct charname_t {
	std::string name;
	float x, y, sx, sy;
	std::vector<imgname_t> imgs;
};

struct image_t {
	std::string name;
	float size;
};

std::map<int, std::string> strtable = { {0, "AABB"},
			{5, "AnimationClip"},
			{19, "AnimationCurve"},
			{34, "AnimationState"},
			{49, "Array"},
			{55, "Base"},
			{60, "BitField"},
			{69, "bitset"},
			{76, "bool"},
			{81, "char"},
			{86, "ColorRGBA"},
			{96, "Component"},
			{106, "data"},
			{111, "deque"},
			{117, "double"},
			{124, "dynamic_array"},
			{138, "FastPropertyName"},
			{155, "first"},
			{161, "float"},
			{167, "Font"},
			{172, "GameObject"},
			{183, "Generic Mono"},
			{196, "GradientNEW"},
			{208, "GUID"},
			{213, "GUIStyle"},
			{222, "int"},
			{226, "list"},
			{231, "long long"},
			{241, "map"},
			{245, "Matrix4x4f"},
			{256, "MdFour"},
			{263, "MonoBehaviour"},
			{277, "MonoScript"},
			{288, "m_ByteSize"},
			{299, "m_Curve"},
			{307, "m_EditorClassIdentifier"},
			{331, "m_EditorHideFlags"},
			{349, "m_Enabled"},
			{359, "m_ExtensionPtr"},
			{374, "m_GameObject"},
			{387, "m_Index"},
			{395, "m_IsArray"},
			{405, "m_IsStatic"},
			{416, "m_MetaFlag"},
			{427, "m_Name"},
			{434, "m_ObjectHideFlags"},
			{452, "m_PrefabInternal"},
			{469, "m_PrefabParentObject"},
			{490, "m_Script"},
			{499, "m_StaticEditorFlags"},
			{519, "m_Type"},
			{526, "m_Version"},
			{536, "Object"},
			{543, "pair"},
			{548, "PPtr<Component>"},
			{564, "PPtr<GameObject>"},
			{581, "PPtr<Material>"},
			{596, "PPtr<MonoBehaviour>"},
			{616, "PPtr<MonoScript>"},
			{633, "PPtr<Object>"},
			{646, "PPtr<Prefab>"},
			{659, "PPtr<Sprite>"},
			{672, "PPtr<TextAsset>"},
			{688, "PPtr<Texture>"},
			{702, "PPtr<Texture2D>"},
			{718, "PPtr<Transform>"},
			{734, "Prefab"},
			{741, "Quaternionf"},
			{753, "Rectf"},
			{759, "RectInt"},
			{767, "RectOffset"},
			{778, "second"},
			{785, "set"},
			{789, "short"},
			{795, "size"},
			{800, "SInt16"},
			{807, "SInt32"},
			{814, "SInt64"},
			{821, "SInt8"},
			{827, "staticvector"},
			{840, "string"},
			{847, "TextAsset"},
			{857, "TextMesh"},
			{866, "Texture"},
			{874, "Texture2D"},
			{884, "Transform"},
			{894, "TypelessData"},
			{907, "UInt16"},
			{914, "UInt32"},
			{921, "UInt64"},
			{928, "UInt8"},
			{934, "unsigned int"},
			{947, "unsigned long long"},
			{966, "unsigned short"},
			{981, "vector"},
			{988, "Vector2f"},
			{997, "Vector3f"},
			{1006, "Vector4f"},
			{1015, "m_ScriptingClassIdentifier"},
			{1042, "Gradient"},
			{1051, "Type*"},
			{1057, "int2_storage"},
			{1070, "int3_storage"},
			{1083, "BoundsInt"},
			{1093, "m_CorrespondingSourceObject"},
			{1121, "m_PrefabInstance"},
			{ 1138, "m_PrefabAsset" } };

int read(const char* data, int length, int& x) {
	x = 0;
	for (int i = 0; i < 4; i++) {
		x = (x << 8) + ((const unsigned char*)data)[i];
	}
	return 4;
}

int read(const char* data, int length, char* x, int n) {
	memcpy(x, data, n);
	return n;
}

int read_(const char* data, int length, int& x) {
	x = 0;
	for (int i = 0; i < 4; i++) {
		x = (x << 8) + ((const unsigned char*)data)[3 - i];
	}
	return 4;
}

int read16(const char* data, int length, int& x) {
	x = 0;
	for (int i = 0; i < 2; i++) {
		x = (x << 8) + ((const unsigned char*)data)[i];
	}
	return 2;
}

int read16_(const char* data, int length, int& x) {
	x = 0;
	for (int i = 0; i < 2; i++) {
		x = (x << 8) + ((const unsigned char*)data)[1 - i];
	}
	return 2;
}

int read8(const char* data, int length, int& x) {
	x = 0;
	for (int i = 0; i < 1; i++) {
		x = (x << 8) + ((const unsigned char*)data)[i];
	}
	return 1;
}

int read(const char* data, int length, long long& x) {
	x = 0;
	for (int i = 0; i < 8; i++) {
		x = (x << 8) + ((const unsigned char*)data)[i];
	}
	return 8;
}

int read_(const char* data, int length, long long& x) {
	x = 0;
	for (int i = 0; i < 8; i++) {
		x = (x << 8) + ((const unsigned char*)data)[7 - i];
	}
	return 8;
}

int read(const char* data, int length, std::string& s) {
	int i;
	for (i = 0; i < length && data[i]; i++);
	s = std::string(data, i);
	return i + 1;
}

int readheader(const char* data, int length, std::vector<block_t>& blocks, std::vector<file_t>& files) {
	int cur = 16;
	int blockcount;
	block_t block;
	cur += read(data + cur, length - cur, blockcount);
	for (int i = 0; i < blockcount; i++) {
		block.headpos = cur;
		cur += read(data + cur, length - cur, block.uncom);
		cur += read(data + cur, length - cur, block.com);
		cur += read16(data + cur, length - cur, block.flag);
		blocks.push_back(block);
	}
	int filecount;
	file_t file;
	cur += read(data + cur, length - cur, filecount);
	for (int i = 0; i < filecount; i++) {
		cur += read(data + cur, length - cur, file.offset);
		cur += read(data + cur, length - cur, file.length);
		cur += read(data + cur, length - cur, file.flag);
		cur += read(data + cur, length - cur, file.name);
		files.push_back(file);
	}
	return 0;
}

int readfile(char* data, int length, std::vector<objectInfo_t>& objects) {
	int MetadataSize, FileSize, Version, DataOffset, EnableTypeTree;
	int cur = 0;
	"char_117_myrrh_2";//17
	"char_117_myrrh_1[alpha]";//24
	"char_117_myrrh_wild#1";//22
	cur += read(data + cur, length - cur, MetadataSize);
	cur += read(data + cur, length - cur, FileSize);
	cur += read(data + cur, length - cur, Version);
	cur += read(data + cur, length - cur, DataOffset);
	int Endianess, TargetPlatform;
	cur += read8(data + cur, length - cur, Endianess);
	cur += 3;
	std::string unityVersion;
	cur += read(data + cur, length - cur, unityVersion);
	cur += read_(data + cur, length - cur, TargetPlatform);
	cur += read8(data + cur, length - cur, EnableTypeTree);
	int typeCount;
	cur += read_(data + cur, length - cur, typeCount);
	std::vector<type_t> types;
	for (int i = 0; i < typeCount; i++) {
		type_t type;
		cur += read_(data + cur, length - cur, type.classID);
		cur += read8(data + cur, length - cur, type.IsStrippedType);
		cur += read16_(data + cur, length - cur, type.ScriptTypeIndex);
		if (type.classID == 114) {
			cur += read(data + cur, length - cur, type.ScriptID, 16);
		}
		cur += read(data + cur, length - cur, type.OldTypeHash, 16);
		if (EnableTypeTree) {
			int numberOfNodes, stringBufferSize;
			cur += read_(data + cur, length - cur, numberOfNodes);
			cur += read_(data + cur, length - cur, stringBufferSize);
			int nodeSize = 24;
			const char* strpos = data + cur + numberOfNodes * nodeSize;
			for (int j = 0; j < numberOfNodes; j++) {
				typeTreeNode_t typeTreeNode;
				cur += read16_(data + cur, length - cur, typeTreeNode.Version);
				cur += read8(data + cur, length - cur, typeTreeNode.Level);
				cur += read8(data + cur, length - cur, typeTreeNode.IsArray);
				cur += read_(data + cur, length - cur, typeTreeNode.TypeStrOffset);
				cur += read_(data + cur, length - cur, typeTreeNode.NameStrOffset);
				cur += read_(data + cur, length - cur, typeTreeNode.ByteSize);
				cur += read_(data + cur, length - cur, typeTreeNode.Index);
				cur += read_(data + cur, length - cur, typeTreeNode.MetaFlag);
				if (typeTreeNode.TypeStrOffset & 0x80000000) {
					typeTreeNode.Type = strtable[typeTreeNode.TypeStrOffset & 0x7fffffff];
				}
				else {
					read(strpos + (typeTreeNode.TypeStrOffset & 0x7fffffff), stringBufferSize - (typeTreeNode.TypeStrOffset & 0x7fffffff), typeTreeNode.Type);
				}
				if (typeTreeNode.NameStrOffset & 0x80000000) {
					typeTreeNode.Name = strtable[typeTreeNode.NameStrOffset & 0x7fffffff];
				}
				else {
					read(strpos + (typeTreeNode.NameStrOffset & 0x7fffffff), stringBufferSize - (typeTreeNode.NameStrOffset & 0x7fffffff), typeTreeNode.Name);
				}
				type.typeTree.push_back(typeTreeNode);
			}
			cur += stringBufferSize;
		}
		types.push_back(type);
	}
	int objectCount;
	cur += read_(data + cur, length - cur, objectCount);
	for (int i = 0; i < objectCount; i++) {
		objectInfo_t objectInfo;
		cur = (cur + 3) / 4 * 4;
		cur += read_(data + cur, length - cur, objectInfo.PathID);
		cur += read_(data + cur, length - cur, objectInfo.byteStart);
		objectInfo.byteStart += DataOffset;
		cur += read_(data + cur, length - cur, objectInfo.byteSize);
		cur += read_(data + cur, length - cur, objectInfo.typeID);
		objectInfo.type = types[objectInfo.typeID];
		objects.push_back(objectInfo);
	}
	return 0;
}

int readtexture2D(const char* data, int length) {
	return 0;
}

int readMonoBehaviour(const char* data, int length, charname_t& charname) {
	int size, cur;
	cur = 32;
	cur += read_(data + cur, length - cur, size);
	int fileid, strl;
	for (int i = 0; i < size; i++) {
		imgname_t imgname;
		cur += read_(data + cur, length - cur, fileid);
		cur += read_(data + cur, length - cur, imgname.rgbid);
		cur += read_(data + cur, length - cur, fileid);
		cur += read_(data + cur, length - cur, imgname.aid);
		cur += read_(data + cur, length - cur, strl);
		imgname.name = std::string(data + cur, strl);
		cur += (strl + 3) / 4 * 4;
		charname.imgs.push_back(imgname);
		if (strl > 0) {
			std::cout << imgname.name << "\n";
		}
	}
	return 0;
}

int change(int x) {
	int y = 0;
	for (int i = 0; i < 4; i++) {
		y = (y << 8) + ((x >> (i * 8)) & 0xff);
	}
	return y;
}

int readts(const char* data, int length, charname_t& ts) {
	memcpy((&ts.x), data + 84, 16);
	return 0;
}

int findobj(std::vector<objectInfo_t>& objects, long long id) {
	for (int k = 0; k < objects.size(); k++) {
		if (objects[k].PathID == id) {
			return k;
		}
	}
	return 0;
}

int dataexp(char* data, std::string out) {

	int sl = ((int*)data)[0];
	FILE* fiout;
	if (!fopen_s(&fiout, (out + "/" + std::string(data + 4, sl) + ".txt").c_str(), "wb")) {
		int fl = ((int*)data)[1 + (sl + 3) / 4];
		fwrite(data + 8 + (sl + 3) / 4 * 4, 1, fl, fiout);
		fclose(fiout);
	}
	return 0;
}

int readfile(char* filename,  std::string dir, std::string out) {
	int i;
	for (i = 0; i < 260; i++) {
		if (filename[i] == '.') {
			break;
		}
	}
	FILE* fi;
	FILE* fiout;
	char* ficon;
	int ficount;
	int cur = 0;
	std::vector<block_t> blocks;
	std::vector<file_t> files;
	if (!fopen_s(&fi, ( dir + "/" + filename).c_str(), "rb")) {
		fseek(fi, 0, SEEK_END);
		ficount = ftell(fi);
		fseek(fi, 0, SEEK_SET);
		ficon = (char*)malloc(ficount * sizeof(char));
		ficount = fread(ficon, 1, ficount, fi);
		std::string s1, s2, s3;
		int format;
		long long bundle;
		int com, uncom, flag;
		cur += read(ficon + cur, ficount - cur, s1);
		cur += read(ficon + cur, ficount - cur, format);
		cur += read(ficon + cur, ficount - cur, s2);
		cur += read(ficon + cur, ficount - cur, s3);
		cur += read(ficon + cur, ficount - cur, bundle);
		cur += read(ficon + cur, ficount - cur, com);
		cur += read(ficon + cur, ficount - cur, uncom);
		cur += read(ficon + cur, ficount - cur, flag);
		char* header;
		header = (char*)malloc(uncom * sizeof(char));
		LZ4_decompress_safe(ficon + cur, header, com, uncom);
		//LZ4_decompress_safe(ficon + cur, header, com, uncom);
		cur += com;
		readheader(header, uncom, blocks, files);
		char* data;
		int datalength = 0;
		int datacur = 0;
		for (i = 0; i < blocks.size(); i++) {
			datalength += blocks[i].uncom;
		}
		data = (char*)malloc(datalength * sizeof(char));
		for (i = 0; i < blocks.size(); i++) {
			LZ4_decompress_safe(ficon + cur, data + datacur, blocks[i].com, blocks[i].uncom);
			cur += blocks[i].com;
			datacur += blocks[i].uncom;
		}
		for (i = 0; i < files.size(); i++) {
			std::vector<objectInfo_t> objects;
			std::cout << files[i].name << "\n";
			readfile(data + files[i].offset, files[i].length, objects);
			for (int j = 0; j < objects.size(); j++) {
				if (objects[j].type.classID == 1) {
					int size = *(int*)(data + files[i].offset + objects[j].byteStart);
					int strlen = ((int*)(data + files[i].offset + objects[j].byteStart))[size * 3 + 2];
					std::string name = std::string(data + files[i].offset + objects[j].byteStart + 4 + size * 12 + 8, strlen);
					if (name == "Front") {
						long long pathid = *(long long*)(data + files[i].offset + objects[j].byteStart + 4 + 3 * 12 + 4);
						int saindex = findobj(objects, pathid);
						long long pathid_SkeletonData = *(long long*)(data + files[i].offset + objects[saindex].byteStart + 9*4);
						int sdindex = findobj(objects, pathid_SkeletonData);

						int strlen = ((int*)(data + files[i].offset + objects[sdindex].byteStart))[7];
						std::string skinname = std::string(data + files[i].offset + objects[sdindex].byteStart + 32, strlen);

						long long pathid_ad = *(long long*)(data + files[i].offset + objects[sdindex].byteStart + 32 + (strlen + 3) / 4 * 4 + 8);
						int adindex = findobj(objects, pathid_ad);
						int astrlen = ((int*)(data + files[i].offset + objects[adindex].byteStart))[7];
						std::string aname = std::string(data + files[i].offset + objects[adindex].byteStart + 32, astrlen);
						long long pathid_at = *(long long*)(data + files[i].offset + objects[adindex].byteStart + 32 + (astrlen + 3) / 4 * 4 + 4);
						int atindex = findobj(objects, pathid_at);

						long long pathid_st = *(long long*)(data + files[i].offset + objects[sdindex].byteStart + 32 + (strlen + 3) / 4 * 4 + 24);
						int stindex = findobj(objects, pathid_st);

						dataexp(data + files[i].offset + objects[stindex].byteStart, out);
						dataexp(data + files[i].offset + objects[atindex].byteStart, out);
						std::cout << skinname << "," << aname << "\n";
					}
				}
			}
		}
	}
	return 0;
}

std::map<std::string, int> namemap;

struct an_t {
	std::string name;
	std::vector<float> dur;
};
std::vector<an_t> ans;

int readan(const char* name,std::string dir) {
	std::string anmis(name, strlen(name) - 10);
	if (anmis.find("#") != std::string::npos) {
		return 0;
	}
	an_t an;
	an.dur = std::vector<float>(namemap.size(), -1.0f);
	an.name = anmis;
	spAtlas* atlas = spAtlas_createFromFile((dir + "/" + anmis + ".atlas.txt").c_str(), 0);
	spSkeletonBinary* binary = spSkeletonBinary_create(atlas);
	spSkeletonData* skeletonData = spSkeletonBinary_readSkeletonDataFile(binary, (dir + "/" + anmis + ".skel.txt").c_str());
	spSkeletonBinary_dispose(binary);
	spAnimationStateData* animationStateData = spAnimationStateData_create(skeletonData);
	spSkeleton* skeleton = spSkeleton_create(skeletonData);
	for (int i = 0; i < skeletonData->animationsCount; i++) {
		std::string name = skeletonData->animations[i]->name;
		if (namemap.find(name) == namemap.end()) {
			namemap.insert(std::make_pair(name, int(namemap.size())));
			an.dur.push_back(0);
		}
		an.dur[namemap[name]] = skeletonData->animations[i]->duration;
	}
	ans.push_back(an);
	return 0;
}

int main(int argc, char** argv)
{
	std::string dir, out;
	if (argc > 1) {
		dir = argv[0];
	}
	else {
		dir = "c:/arknights-hg-1001/assets/AB/Android/charpack";
	}
	if (argc > 2) {
		out = argv[1];
	}
	else {
		out = "c:/spine";
	}
	intptr_t handle;
	struct _finddata_t FileInfo;
	handle = _findfirst((dir + "/*.ab").c_str(), &FileInfo);
	if (handle != -1) {
		printf("%s\n", FileInfo.name);
		readfile(FileInfo.name, dir,out);
		while (!_findnext(handle, &FileInfo)) {
			printf("%s\n", FileInfo.name);
			readfile(FileInfo.name, dir,out);
		}
		_findclose(handle);
	}
  
	if (argc > 2) {
		dir = argv[1];
	}
	else {
		dir = "c:/spine";
	}
	if (argc > 3) {
		out = argv[2];
	}
	else {
		out = "c:/animf.html";
	}
	FILE* fiout;
	if (!fopen_s(&fiout, out.c_str(), "wb")) {
		intptr_t handle;
		struct _finddata_t FileInfo;
		handle = _findfirst((dir + "/*.atlas.txt").c_str(), &FileInfo);
		if (handle != -1) {
			printf("%s\n", FileInfo.name);
			readan(FileInfo.name, dir);
			while (!_findnext(handle, &FileInfo)) {
				printf("%s\n", FileInfo.name);
				readan(FileInfo.name, dir);
			}
			_findclose(handle);
		}
		std::vector<std::string> names(namemap.size());
		for (auto it = namemap.begin(); it != namemap.end(); it++) {
			names[it->second] = it->first;
		}

		fprintf(fiout, "<html><head><meta charset=\"UTF-8\"></head><body><table rules=\"all\"><tr>");
		fprintf(fiout, "<td></td>");
		for (int i = 0; i < names.size(); i++) {
			fprintf(fiout, "<td>%s</td>", names[i].c_str());
		}
		fprintf(fiout, "</tr>");
		for (int i = 0; i < ans.size(); i++) {

			fprintf(fiout, "<tr>");
			fprintf(fiout, "<td>%s</td>", ans[i].name.c_str());
			int j = 0;
			for (; j < ans[i].dur.size(); j++) {
				if (ans[i].dur[j] != -1) {
					fprintf(fiout, "<td>%f</td>", ans[i].dur[j]);
				}
				else {
					fprintf(fiout, "<td></td>");
				}
			}
			for (; j < namemap.size(); j++) {
				fprintf(fiout, "<td></td>");
			}
			fprintf(fiout, "</tr>");
		}
		fprintf(fiout, "</table></body></html>");
		fclose(fiout);
	}
}


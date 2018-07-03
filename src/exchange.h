void createexportlist(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,grid_t grid,struct Zoltan_Struct *ztn,celltype_t* celltype,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata);
void exchangecleanup(systeminfo_t systeminfo,cellsinfo_t cellsinfo,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata);
void cellssendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void cellswait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datasendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datawait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);

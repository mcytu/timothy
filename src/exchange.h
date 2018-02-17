void createexportlist(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,celltype_t* celltype,cellcommdata_t *cellcommdata);
void exchangecleanup(systeminfo_t systeminfo,cellsinfo_t cellsinfo,cellcommdata_t *cellcommdata);
void cellssendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void cellswait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datasendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datawait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);

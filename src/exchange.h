void createexportlist(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,celltype_t* celltype,commdata_t *commdata);
void exchangecleanup(systeminfo_t systeminfo,cellsinfo_t cellsinfo,commdata_t *commdata);
void cellssendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, commdata_t *commdata);
void cellswait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, commdata_t *commdata);
void datasendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, commdata_t *commdata);
void datawait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, commdata_t *commdata);

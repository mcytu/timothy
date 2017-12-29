void createexportlist(system_t system,settings_t settings,cellsinfo_t cellsinfo,celltype_t* celltype,commdata_t *commdata);
void exchangecleanup(system_t system,cellsinfo_t cellsinfo,commdata_t *commdata);
void cellssendrecv(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata);
void cellswait(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata);
void datasendrecv(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata);
void datawait(system_t system, cellsinfo_t cellsinfo, commdata_t *commdata);

void getsysteminfo(systeminfo_t* systeminfo);
void initialsettings(settings_t* settings);
void initialcelltype(int numberofcelltypes,int numberoffields,celltype_t* celltype);
void initialfields(int numberoffields,environment_t* environment);
void initialisation(int argc, char **argv, systeminfo_t *systeminfo, settings_t* settings,celltype_t** celltype,environment_t** environment);
void initcount(cellcount_t *cellcount);
void allocatecells(systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo);
void printinfo(systeminfo_t systeminfo);

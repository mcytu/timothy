void getsysteminfo(system_t* system);
void initialsettings(settings_t* settings);
void initialcelltype(int numberofcelltypes,int numberoffields,celltype_t* celltype);
void initialfields(int numberoffields,environment_t* environment);
void initialisation(int argc, char **argv, system_t *system, settings_t* settings,celltype_t** celltype,environment_t** environment);
void initcount(cellcount_t *cellcount);
void allocatecells(system_t system,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo);
void printinfo(system_t system);

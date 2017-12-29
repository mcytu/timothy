int checkendiannes(system_t *system);
void swapnbyte(char *data, int n, int m);
void updateglobalcounts(cellsinfo_t* cellsinfo);
void terminate(system_t system, char *msg, char *file, int line);
void stopRun(int ierr, char *name, char *file, int line);
size_t getMemoryPerProcess(int32_t lsize);
void getLocalRankAndSize(int rank, int size, int32_t * lrank, int32_t * lsize);
void randomstreaminit(system_t *system,settings_t *settings);
void printstatistics(system_t system,settings_t settings,cellsinfo_t cellsinfo,statistics_t* statistics);

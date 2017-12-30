int checkendiannes(systeminfo_t *systeminfo);
void swapnbyte(char *data, int n, int m);
void updateglobalcounts(cellsinfo_t* cellsinfo);
void terminate(systeminfo_t systeminfo, char *msg, char *file, int line);
void stopRun(int ierr, char *name, char *file, int line);
size_t getMemoryPerProcess(int32_t lsize);
void getLocalRankAndSize(int rank, int size, int32_t * lrank, int32_t * lsize);
void randomstreaminit(systeminfo_t *systeminfo,settings_t *settings);
void printstatistics(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,statistics_t* statistics);

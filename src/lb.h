void lbinit(int argc, char **argv, MPI_Comm Comm,systeminfo_t systeminfo,struct Zoltan_Struct **ztn,cellsinfo_t *cellsinfo);
void lbexchange(systeminfo_t systeminfo,struct Zoltan_Struct *ztn);
void lbdestroy(struct Zoltan_Struct **ztn);

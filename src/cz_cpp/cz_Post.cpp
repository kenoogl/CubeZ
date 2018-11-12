/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

#include "cz.h"

int CZ::Post()
{
  if (!fph) fclose(fph);


  char tmp_fname[30];
  int gc = GUIDE;

  if (debug_mode==1) {
    sprintf( tmp_fname, "p_%05d.sph", myRank );
    fileout_(size, &gc, P, pitch, origin, tmp_fname);

    exact_(size, &gc, ERR, pitch, origin);
    double errmax = 0.0;
    err_  (size, &gc, &errmax, P, ERR);
    sprintf( tmp_fname, "e_%05d.sph", myRank );
    fileout_(size, &gc, ERR, pitch, origin, tmp_fname);
  }


  FILE *fp = NULL;

  Hostonly_
  {
    if ( !(fp=fopen("profiling.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
      Exit(0);
    }
  }


  // 測定結果の集計(gatherメソッドは全ノードで呼ぶこと)
  PM.gather();

  char str[100];
  sprintf(str, "CubeZ %s", CZ_VERSION);

  string HostName = GetHostName();

  // 結果出力(排他測定のみ)
  PM.print(stdout, HostName, str);
  PM.print(fp, HostName, str);
  PM.printDetail(fp);

  for (int i=0; i<numProc; i++)
  {
    PM.printThreads(fp, i, 0); // time cost order
  }

	PM.printLegend(fp);

  Hostonly_ {
    fflush(fp);
    fclose(fp);
  }

  return 1;

}

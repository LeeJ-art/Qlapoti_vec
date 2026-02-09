#include <stdint.h>
#include <stdio.h>
#include <arm_neon.h>
#include <hd.h>
#include <string.h>

void initial_montback_array(fp_t* mb, fp_t* imb){
    //2^522
    mb[0][0] = 155344;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 0;
    mb[0][7] = 0;
    mb[0][8] = 281474976710656;
    //2^(513*2-522)
    imb[0][0] = 0;
    imb[0][1] = 0;
    imb[0][2] = 0;
    imb[0][3] = 0;
    imb[0][4] = 0;
    imb[0][5] = 0;
    imb[0][6] = 0;
    imb[0][7] = 0;
    imb[0][8] = 281474976710656;
}

void initial_q_value_1(uint32_t* q1){
    uint32_t tmp[FP_LIMBS] = {155344, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2048};
    memmove(q1, tmp, sizeof(uint32_t)*FP_LIMBS);
}

void ec_montback_array(fp_t* mb){
    //2^(522*4 - 513*3)
    mb[0][0] = 20849998274863;
    mb[0][1] = 0;
    mb[0][2] = 0;
    mb[0][3] = 0;
    mb[0][4] = 0;
    mb[0][5] = 0;
    mb[0][6] = 0;
    mb[0][7] = 0;
    mb[0][8] = 193514046488576;
}


void transpose(uint32x4_t *Out, theta_point_t In){
    uint64_t mask = ((uint64_t)1<<PER_LIMB)-1;
    uint32_t in32[4][FP2_LIMBS] = {0};

    // re-cast
    in32[0][0] = (uint32_t)(mask & In.x.re[0]);
    in32[1][0] = (uint32_t)(mask & In.y.re[0]);
    in32[2][0] = (uint32_t)(mask & In.z.re[0]);
    in32[3][0] = (uint32_t)(mask & In.t.re[0]);

    in32[0][1] = (uint32_t)(mask & ((In.x.re[0]>>29) + (In.x.re[1]<<28)));
    in32[1][1] = (uint32_t)(mask & ((In.y.re[0]>>29) + (In.y.re[1]<<28)));
    in32[2][1] = (uint32_t)(mask & ((In.z.re[0]>>29) + (In.z.re[1]<<28)));
    in32[3][1] = (uint32_t)(mask & ((In.t.re[0]>>29) + (In.t.re[1]<<28)));

    in32[0][2] = (uint32_t)(mask & (In.x.re[1]>>1));
    in32[1][2] = (uint32_t)(mask & (In.y.re[1]>>1));
    in32[2][2] = (uint32_t)(mask & (In.z.re[1]>>1));
    in32[3][2] = (uint32_t)(mask & (In.t.re[1]>>1));

    in32[0][3] = (uint32_t)(mask & ((In.x.re[1]>>30) + (In.x.re[2]<<27)));
    in32[1][3] = (uint32_t)(mask & ((In.y.re[1]>>30) + (In.y.re[2]<<27)));
    in32[2][3] = (uint32_t)(mask & ((In.z.re[1]>>30) + (In.z.re[2]<<27)));
    in32[3][3] = (uint32_t)(mask & ((In.t.re[1]>>30) + (In.t.re[2]<<27)));

    in32[0][4] = (uint32_t)(mask & (In.x.re[2]>>2));
    in32[1][4] = (uint32_t)(mask & (In.y.re[2]>>2));
    in32[2][4] = (uint32_t)(mask & (In.z.re[2]>>2));
    in32[3][4] = (uint32_t)(mask & (In.t.re[2]>>2));

    in32[0][5] = (uint32_t)(mask & ((In.x.re[2]>>31) + (In.x.re[3]<<26)));
    in32[1][5] = (uint32_t)(mask & ((In.y.re[2]>>31) + (In.y.re[3]<<26)));
    in32[2][5] = (uint32_t)(mask & ((In.z.re[2]>>31) + (In.z.re[3]<<26)));
    in32[3][5] = (uint32_t)(mask & ((In.t.re[2]>>31) + (In.t.re[3]<<26)));

    in32[0][6] = (uint32_t)(mask & (In.x.re[3]>>3));
    in32[1][6] = (uint32_t)(mask & (In.y.re[3]>>3));
    in32[2][6] = (uint32_t)(mask & (In.z.re[3]>>3));
    in32[3][6] = (uint32_t)(mask & (In.t.re[3]>>3));

    in32[0][7] = (uint32_t)(mask & ((In.x.re[3]>>32) + (In.x.re[4]<<25)));
    in32[1][7] = (uint32_t)(mask & ((In.y.re[3]>>32) + (In.y.re[4]<<25)));
    in32[2][7] = (uint32_t)(mask & ((In.z.re[3]>>32) + (In.z.re[4]<<25)));
    in32[3][7] = (uint32_t)(mask & ((In.t.re[3]>>32) + (In.t.re[4]<<25)));

    in32[0][8] = (uint32_t)(mask & (In.x.re[4]>>4));
    in32[1][8] = (uint32_t)(mask & (In.y.re[4]>>4));
    in32[2][8] = (uint32_t)(mask & (In.z.re[4]>>4));
    in32[3][8] = (uint32_t)(mask & (In.t.re[4]>>4));

    in32[0][9] = (uint32_t)(mask & ((In.x.re[4]>>33) + (In.x.re[5]<<24)));
    in32[1][9] = (uint32_t)(mask & ((In.y.re[4]>>33) + (In.y.re[5]<<24)));
    in32[2][9] = (uint32_t)(mask & ((In.z.re[4]>>33) + (In.z.re[5]<<24)));
    in32[3][9] = (uint32_t)(mask & ((In.t.re[4]>>33) + (In.t.re[5]<<24)));

    in32[0][10] = (uint32_t)(mask & (In.x.re[5]>>5));
    in32[1][10] = (uint32_t)(mask & (In.y.re[5]>>5));
    in32[2][10] = (uint32_t)(mask & (In.z.re[5]>>5));
    in32[3][10] = (uint32_t)(mask & (In.t.re[5]>>5));

    in32[0][11] = (uint32_t)(mask & ((In.x.re[5]>>34) + (In.x.re[6]<<23)));
    in32[1][11] = (uint32_t)(mask & ((In.y.re[5]>>34) + (In.y.re[6]<<23)));
    in32[2][11] = (uint32_t)(mask & ((In.z.re[5]>>34) + (In.z.re[6]<<23)));
    in32[3][11] = (uint32_t)(mask & ((In.t.re[5]>>34) + (In.t.re[6]<<23)));

    in32[0][12] = (uint32_t)(mask & (In.x.re[6]>>6));
    in32[1][12] = (uint32_t)(mask & (In.y.re[6]>>6));
    in32[2][12] = (uint32_t)(mask & (In.z.re[6]>>6));
    in32[3][12] = (uint32_t)(mask & (In.t.re[6]>>6));
    
    in32[0][13] = (uint32_t)(mask & ((In.x.re[6]>>35) + (In.x.re[7]<<22)));
    in32[1][13] = (uint32_t)(mask & ((In.y.re[6]>>35) + (In.y.re[7]<<22)));
    in32[2][13] = (uint32_t)(mask & ((In.z.re[6]>>35) + (In.z.re[7]<<22)));
    in32[3][13] = (uint32_t)(mask & ((In.t.re[6]>>35) + (In.t.re[7]<<22)));

    in32[0][14] = (uint32_t)(mask & (In.x.re[7]>>7));
    in32[1][14] = (uint32_t)(mask & (In.y.re[7]>>7));
    in32[2][14] = (uint32_t)(mask & (In.z.re[7]>>7));
    in32[3][14] = (uint32_t)(mask & (In.t.re[7]>>7));
    
    in32[0][15] = (uint32_t)(mask & ((In.x.re[7]>>36) + (In.x.re[8]<<21)));
    in32[1][15] = (uint32_t)(mask & ((In.y.re[7]>>36) + (In.y.re[8]<<21)));
    in32[2][15] = (uint32_t)(mask & ((In.z.re[7]>>36) + (In.z.re[8]<<21)));
    in32[3][15] = (uint32_t)(mask & ((In.t.re[7]>>36) + (In.t.re[8]<<21)));

    in32[0][16] = (uint32_t)(mask & (In.x.re[8]>>8));
    in32[1][16] = (uint32_t)(mask & (In.y.re[8]>>8));
    in32[2][16] = (uint32_t)(mask & (In.z.re[8]>>8));
    in32[3][16] = (uint32_t)(mask & (In.t.re[8]>>8));
    
    in32[0][17] = (uint32_t)(mask & (In.x.re[8]>>37));
    in32[1][17] = (uint32_t)(mask & (In.y.re[8]>>37));
    in32[2][17] = (uint32_t)(mask & (In.z.re[8]>>37));
    in32[3][17] = (uint32_t)(mask & (In.t.re[8]>>37));

    // im-cast
    in32[0][18] = (uint32_t)(mask & In.x.im[0]);
    in32[1][18] = (uint32_t)(mask & In.y.im[0]);
    in32[2][18] = (uint32_t)(mask & In.z.im[0]);
    in32[3][18] = (uint32_t)(mask & In.t.im[0]);

    in32[0][19] = (uint32_t)(mask & ((In.x.im[0]>>29) + (In.x.im[1]<<28)));
    in32[1][19] = (uint32_t)(mask & ((In.y.im[0]>>29) + (In.y.im[1]<<28)));
    in32[2][19] = (uint32_t)(mask & ((In.z.im[0]>>29) + (In.z.im[1]<<28)));
    in32[3][19] = (uint32_t)(mask & ((In.t.im[0]>>29) + (In.t.im[1]<<28)));

    in32[0][20] = (uint32_t)(mask & (In.x.im[1]>>1));
    in32[1][20] = (uint32_t)(mask & (In.y.im[1]>>1));
    in32[2][20] = (uint32_t)(mask & (In.z.im[1]>>1));
    in32[3][20] = (uint32_t)(mask & (In.t.im[1]>>1));

    in32[0][21] = (uint32_t)(mask & ((In.x.im[1]>>30) + (In.x.im[2]<<27)));
    in32[1][21] = (uint32_t)(mask & ((In.y.im[1]>>30) + (In.y.im[2]<<27)));
    in32[2][21] = (uint32_t)(mask & ((In.z.im[1]>>30) + (In.z.im[2]<<27)));
    in32[3][21] = (uint32_t)(mask & ((In.t.im[1]>>30) + (In.t.im[2]<<27)));

    in32[0][22] = (uint32_t)(mask & (In.x.im[2]>>2));
    in32[1][22] = (uint32_t)(mask & (In.y.im[2]>>2));
    in32[2][22] = (uint32_t)(mask & (In.z.im[2]>>2));
    in32[3][22] = (uint32_t)(mask & (In.t.im[2]>>2));

    in32[0][23] = (uint32_t)(mask & ((In.x.im[2]>>31) + (In.x.im[3]<<26)));
    in32[1][23] = (uint32_t)(mask & ((In.y.im[2]>>31) + (In.y.im[3]<<26)));
    in32[2][23] = (uint32_t)(mask & ((In.z.im[2]>>31) + (In.z.im[3]<<26)));
    in32[3][23] = (uint32_t)(mask & ((In.t.im[2]>>31) + (In.t.im[3]<<26)));

    in32[0][24] = (uint32_t)(mask & (In.x.im[3]>>3));
    in32[1][24] = (uint32_t)(mask & (In.y.im[3]>>3));
    in32[2][24] = (uint32_t)(mask & (In.z.im[3]>>3));
    in32[3][24] = (uint32_t)(mask & (In.t.im[3]>>3));

    in32[0][25] = (uint32_t)(mask & ((In.x.im[3]>>32) + (In.x.im[4]<<25)));
    in32[1][25] = (uint32_t)(mask & ((In.y.im[3]>>32) + (In.y.im[4]<<25)));
    in32[2][25] = (uint32_t)(mask & ((In.z.im[3]>>32) + (In.z.im[4]<<25)));
    in32[3][25] = (uint32_t)(mask & ((In.t.im[3]>>32) + (In.t.im[4]<<25)));

    in32[0][26] = (uint32_t)(mask & (In.x.im[4]>>4));
    in32[1][26] = (uint32_t)(mask & (In.y.im[4]>>4));
    in32[2][26] = (uint32_t)(mask & (In.z.im[4]>>4));
    in32[3][26] = (uint32_t)(mask & (In.t.im[4]>>4));

    in32[0][27] = (uint32_t)(mask & ((In.x.im[4]>>33) + (In.x.im[5]<<24)));
    in32[1][27] = (uint32_t)(mask & ((In.y.im[4]>>33) + (In.y.im[5]<<24)));
    in32[2][27] = (uint32_t)(mask & ((In.z.im[4]>>33) + (In.z.im[5]<<24)));
    in32[3][27] = (uint32_t)(mask & ((In.t.im[4]>>33) + (In.t.im[5]<<24)));

    in32[0][28] = (uint32_t)(mask & (In.x.im[5]>>5));
    in32[1][28] = (uint32_t)(mask & (In.y.im[5]>>5));
    in32[2][28] = (uint32_t)(mask & (In.z.im[5]>>5));
    in32[3][28] = (uint32_t)(mask & (In.t.im[5]>>5));

    in32[0][29] = (uint32_t)(mask & ((In.x.im[5]>>34) + (In.x.im[6]<<23)));
    in32[1][29] = (uint32_t)(mask & ((In.y.im[5]>>34) + (In.y.im[6]<<23)));
    in32[2][29] = (uint32_t)(mask & ((In.z.im[5]>>34) + (In.z.im[6]<<23)));
    in32[3][29] = (uint32_t)(mask & ((In.t.im[5]>>34) + (In.t.im[6]<<23)));

    in32[0][30] = (uint32_t)(mask & (In.x.im[6]>>6));
    in32[1][30] = (uint32_t)(mask & (In.y.im[6]>>6));
    in32[2][30] = (uint32_t)(mask & (In.z.im[6]>>6));
    in32[3][30] = (uint32_t)(mask & (In.t.im[6]>>6));
    
    in32[0][31] = (uint32_t)(mask & ((In.x.im[6]>>35) + (In.x.im[7]<<22)));
    in32[1][31] = (uint32_t)(mask & ((In.y.im[6]>>35) + (In.y.im[7]<<22)));
    in32[2][31] = (uint32_t)(mask & ((In.z.im[6]>>35) + (In.z.im[7]<<22)));
    in32[3][31] = (uint32_t)(mask & ((In.t.im[6]>>35) + (In.t.im[7]<<22)));

    in32[0][32] = (uint32_t)(mask & (In.x.im[7]>>7));
    in32[1][32] = (uint32_t)(mask & (In.y.im[7]>>7));
    in32[2][32] = (uint32_t)(mask & (In.z.im[7]>>7));
    in32[3][32] = (uint32_t)(mask & (In.t.im[7]>>7));
    
    in32[0][33] = (uint32_t)(mask & ((In.x.im[7]>>36) + (In.x.im[8]<<21)));
    in32[1][33] = (uint32_t)(mask & ((In.y.im[7]>>36) + (In.y.im[8]<<21)));
    in32[2][33] = (uint32_t)(mask & ((In.z.im[7]>>36) + (In.z.im[8]<<21)));
    in32[3][33] = (uint32_t)(mask & ((In.t.im[7]>>36) + (In.t.im[8]<<21)));

    in32[0][34] = (uint32_t)(mask & (In.x.im[8]>>8));
    in32[1][34] = (uint32_t)(mask & (In.y.im[8]>>8));
    in32[2][34] = (uint32_t)(mask & (In.z.im[8]>>8));
    in32[3][34] = (uint32_t)(mask & (In.t.im[8]>>8));
    
    in32[0][35] = (uint32_t)(mask & (In.x.im[8]>>37));
    in32[1][35] = (uint32_t)(mask & (In.y.im[8]>>37));
    in32[2][35] = (uint32_t)(mask & (In.z.im[8]>>37));
    in32[3][35] = (uint32_t)(mask & (In.t.im[8]>>37));

    // transpose
    for(int i = 0;i<FP2_LIMBS; i++){
        uint32x4_t tmp = {in32[0][i], in32[1][i], in32[2][i], in32[3][i]};
        Out[i] = tmp;
    }
}

void itranspose(theta_point_t *Out, uint32x4_t *In){
    uint32_t in32[4][FP2_LIMBS];

    // itranspose
    for(int i = 0;i<4;i++){
        for(int j = 0;j<FP2_LIMBS;j++){
            in32[i][j] = In[j][i];
        }
    }

    // re-cast
    uint64_t mask = ((uint64_t)1<<57)-1;
    Out->x.re[0] = mask & (((uint64_t)in32[0][0])  + ((uint64_t)in32[0][1] <<29));
    Out->x.im[0] = mask & (((uint64_t)in32[0][18]) + ((uint64_t)in32[0][19]<<29));
    Out->y.re[0] = mask & (((uint64_t)in32[1][0])  + ((uint64_t)in32[1][1] <<29));
    Out->y.im[0] = mask & (((uint64_t)in32[1][18]) + ((uint64_t)in32[1][19]<<29));
    Out->z.re[0] = mask & (((uint64_t)in32[2][0])  + ((uint64_t)in32[2][1] <<29));
    Out->z.im[0] = mask & (((uint64_t)in32[2][18]) + ((uint64_t)in32[2][19]<<29));
    Out->t.re[0] = mask & (((uint64_t)in32[3][0])  + ((uint64_t)in32[3][1] <<29));
    Out->t.im[0] = mask & (((uint64_t)in32[3][18]) + ((uint64_t)in32[3][19]<<29));

    Out->x.re[1] = mask & (((uint64_t)in32[0][1] >>28) + ((uint64_t)in32[0][2] <<1) + ((uint64_t)in32[0][3] <<30));
    Out->x.im[1] = mask & (((uint64_t)in32[0][19]>>28) + ((uint64_t)in32[0][20]<<1) + ((uint64_t)in32[0][21]<<30));
    Out->y.re[1] = mask & (((uint64_t)in32[1][1] >>28) + ((uint64_t)in32[1][2] <<1) + ((uint64_t)in32[1][3] <<30));
    Out->y.im[1] = mask & (((uint64_t)in32[1][19]>>28) + ((uint64_t)in32[1][20]<<1) + ((uint64_t)in32[1][21]<<30));
    Out->z.re[1] = mask & (((uint64_t)in32[2][1] >>28) + ((uint64_t)in32[2][2] <<1) + ((uint64_t)in32[2][3] <<30));
    Out->z.im[1] = mask & (((uint64_t)in32[2][19]>>28) + ((uint64_t)in32[2][20]<<1) + ((uint64_t)in32[2][21]<<30));
    Out->t.re[1] = mask & (((uint64_t)in32[3][1] >>28) + ((uint64_t)in32[3][2] <<1) + ((uint64_t)in32[3][3] <<30));
    Out->t.im[1] = mask & (((uint64_t)in32[3][19]>>28) + ((uint64_t)in32[3][20]<<1) + ((uint64_t)in32[3][21]<<30));

    Out->x.re[2] = mask & (((uint64_t)in32[0][3] >>27) + ((uint64_t)in32[0][4] <<2) + ((uint64_t)in32[0][5] <<31));
    Out->x.im[2] = mask & (((uint64_t)in32[0][21]>>27) + ((uint64_t)in32[0][22]<<2) + ((uint64_t)in32[0][23]<<31));
    Out->y.re[2] = mask & (((uint64_t)in32[1][3] >>27) + ((uint64_t)in32[1][4] <<2) + ((uint64_t)in32[1][5] <<31));
    Out->y.im[2] = mask & (((uint64_t)in32[1][21]>>27) + ((uint64_t)in32[1][22]<<2) + ((uint64_t)in32[1][23]<<31));
    Out->z.re[2] = mask & (((uint64_t)in32[2][3] >>27) + ((uint64_t)in32[2][4] <<2) + ((uint64_t)in32[2][5] <<31));
    Out->z.im[2] = mask & (((uint64_t)in32[2][21]>>27) + ((uint64_t)in32[2][22]<<2) + ((uint64_t)in32[2][23]<<31));
    Out->t.re[2] = mask & (((uint64_t)in32[3][3] >>27) + ((uint64_t)in32[3][4] <<2) + ((uint64_t)in32[3][5] <<31));
    Out->t.im[2] = mask & (((uint64_t)in32[3][21]>>27) + ((uint64_t)in32[3][22]<<2) + ((uint64_t)in32[3][23]<<31));

    Out->x.re[3] = mask & (((uint64_t)in32[0][5] >>26) + ((uint64_t)in32[0][6] <<3) + ((uint64_t)in32[0][7] <<32));
    Out->x.im[3] = mask & (((uint64_t)in32[0][23]>>26) + ((uint64_t)in32[0][24]<<3) + ((uint64_t)in32[0][25]<<32));
    Out->y.re[3] = mask & (((uint64_t)in32[1][5] >>26) + ((uint64_t)in32[1][6] <<3) + ((uint64_t)in32[1][7] <<32));
    Out->y.im[3] = mask & (((uint64_t)in32[1][23]>>26) + ((uint64_t)in32[1][24]<<3) + ((uint64_t)in32[1][25]<<32));
    Out->z.re[3] = mask & (((uint64_t)in32[2][5] >>26) + ((uint64_t)in32[2][6] <<3) + ((uint64_t)in32[2][7] <<32));
    Out->z.im[3] = mask & (((uint64_t)in32[2][23]>>26) + ((uint64_t)in32[2][24]<<3) + ((uint64_t)in32[2][25]<<32));
    Out->t.re[3] = mask & (((uint64_t)in32[3][5] >>26) + ((uint64_t)in32[3][6] <<3) + ((uint64_t)in32[3][7] <<32));
    Out->t.im[3] = mask & (((uint64_t)in32[3][23]>>26) + ((uint64_t)in32[3][24]<<3) + ((uint64_t)in32[3][25]<<32));

    Out->x.re[4] = mask & (((uint64_t)in32[0][7] >>25) + ((uint64_t)in32[0][8] <<4) + ((uint64_t)in32[0][9] <<33));
    Out->x.im[4] = mask & (((uint64_t)in32[0][25]>>25) + ((uint64_t)in32[0][26]<<4) + ((uint64_t)in32[0][27]<<33));
    Out->y.re[4] = mask & (((uint64_t)in32[1][7] >>25) + ((uint64_t)in32[1][8] <<4) + ((uint64_t)in32[1][9] <<33));
    Out->y.im[4] = mask & (((uint64_t)in32[1][25]>>25) + ((uint64_t)in32[1][26]<<4) + ((uint64_t)in32[1][27]<<33));
    Out->z.re[4] = mask & (((uint64_t)in32[2][7] >>25) + ((uint64_t)in32[2][8] <<4) + ((uint64_t)in32[2][9] <<33));
    Out->z.im[4] = mask & (((uint64_t)in32[2][25]>>25) + ((uint64_t)in32[2][26]<<4) + ((uint64_t)in32[2][27]<<33));
    Out->t.re[4] = mask & (((uint64_t)in32[3][7] >>25) + ((uint64_t)in32[3][8] <<4) + ((uint64_t)in32[3][9] <<33));
    Out->t.im[4] = mask & (((uint64_t)in32[3][25]>>25) + ((uint64_t)in32[3][26]<<4) + ((uint64_t)in32[3][27]<<33));

    Out->x.re[5] = mask & (((uint64_t)in32[0][9] >>24) + ((uint64_t)in32[0][10]<<5) + ((uint64_t)in32[0][11] <<34));
    Out->x.im[5] = mask & (((uint64_t)in32[0][27]>>24) + ((uint64_t)in32[0][28]<<5) + ((uint64_t)in32[0][29] <<34));
    Out->y.re[5] = mask & (((uint64_t)in32[1][9] >>24) + ((uint64_t)in32[1][10]<<5) + ((uint64_t)in32[1][11] <<34));
    Out->y.im[5] = mask & (((uint64_t)in32[1][27]>>24) + ((uint64_t)in32[1][28]<<5) + ((uint64_t)in32[1][29] <<34));
    Out->z.re[5] = mask & (((uint64_t)in32[2][9] >>24) + ((uint64_t)in32[2][10]<<5) + ((uint64_t)in32[2][11] <<34));
    Out->z.im[5] = mask & (((uint64_t)in32[2][27]>>24) + ((uint64_t)in32[2][28]<<5) + ((uint64_t)in32[2][29] <<34));
    Out->t.re[5] = mask & (((uint64_t)in32[3][9] >>24) + ((uint64_t)in32[3][10]<<5) + ((uint64_t)in32[3][11] <<34));
    Out->t.im[5] = mask & (((uint64_t)in32[3][27]>>24) + ((uint64_t)in32[3][28]<<5) + ((uint64_t)in32[3][29] <<34));

    Out->x.re[6] = mask & (((uint64_t)in32[0][11] >>23) + ((uint64_t)in32[0][12] <<6) + ((uint64_t)in32[0][13] <<35));
    Out->x.im[6] = mask & (((uint64_t)in32[0][29] >>23) + ((uint64_t)in32[0][30] <<6) + ((uint64_t)in32[0][31] <<35));
    Out->y.re[6] = mask & (((uint64_t)in32[1][11] >>23) + ((uint64_t)in32[1][12] <<6) + ((uint64_t)in32[1][13] <<35));
    Out->y.im[6] = mask & (((uint64_t)in32[1][29] >>23) + ((uint64_t)in32[1][30] <<6) + ((uint64_t)in32[1][31] <<35));
    Out->z.re[6] = mask & (((uint64_t)in32[2][11] >>23) + ((uint64_t)in32[2][12] <<6) + ((uint64_t)in32[2][13] <<35));
    Out->z.im[6] = mask & (((uint64_t)in32[2][29] >>23) + ((uint64_t)in32[2][30] <<6) + ((uint64_t)in32[2][31] <<35));
    Out->t.re[6] = mask & (((uint64_t)in32[3][11] >>23) + ((uint64_t)in32[3][12] <<6) + ((uint64_t)in32[3][13] <<35));
    Out->t.im[6] = mask & (((uint64_t)in32[3][29] >>23) + ((uint64_t)in32[3][30] <<6) + ((uint64_t)in32[3][31] <<35));

    Out->x.re[7] = mask & (((uint64_t)in32[0][13] >>22) + ((uint64_t)in32[0][14] <<7) + ((uint64_t)in32[0][15] <<36));
    Out->x.im[7] = mask & (((uint64_t)in32[0][31] >>22) + ((uint64_t)in32[0][32] <<7) + ((uint64_t)in32[0][33] <<36));
    Out->y.re[7] = mask & (((uint64_t)in32[1][13] >>22) + ((uint64_t)in32[1][14] <<7) + ((uint64_t)in32[1][15] <<36));
    Out->y.im[7] = mask & (((uint64_t)in32[1][31] >>22) + ((uint64_t)in32[1][32] <<7) + ((uint64_t)in32[1][33] <<36));
    Out->z.re[7] = mask & (((uint64_t)in32[2][13] >>22) + ((uint64_t)in32[2][14] <<7) + ((uint64_t)in32[2][15] <<36));
    Out->z.im[7] = mask & (((uint64_t)in32[2][31] >>22) + ((uint64_t)in32[2][32] <<7) + ((uint64_t)in32[2][33] <<36));
    Out->t.re[7] = mask & (((uint64_t)in32[3][13] >>22) + ((uint64_t)in32[3][14] <<7) + ((uint64_t)in32[3][15] <<36));
    Out->t.im[7] = mask & (((uint64_t)in32[3][31] >>22) + ((uint64_t)in32[3][32] <<7) + ((uint64_t)in32[3][33] <<36));

    Out->x.re[8] = mask & (((uint64_t)in32[0][15] >>21) + ((uint64_t)in32[0][16] <<8) + ((uint64_t)in32[0][17] <<37));
    Out->x.im[8] = mask & (((uint64_t)in32[0][33] >>21) + ((uint64_t)in32[0][34] <<8) + ((uint64_t)in32[0][35] <<37));
    Out->y.re[8] = mask & (((uint64_t)in32[1][15] >>21) + ((uint64_t)in32[1][16] <<8) + ((uint64_t)in32[1][17] <<37));
    Out->y.im[8] = mask & (((uint64_t)in32[1][33] >>21) + ((uint64_t)in32[1][34] <<8) + ((uint64_t)in32[1][35] <<37));
    Out->z.re[8] = mask & (((uint64_t)in32[2][15] >>21) + ((uint64_t)in32[2][16] <<8) + ((uint64_t)in32[2][17] <<37));
    Out->z.im[8] = mask & (((uint64_t)in32[2][33] >>21) + ((uint64_t)in32[2][34] <<8) + ((uint64_t)in32[2][35] <<37));
    Out->t.re[8] = mask & (((uint64_t)in32[3][15] >>21) + ((uint64_t)in32[3][16] <<8) + ((uint64_t)in32[3][17] <<37));
    Out->t.im[8] = mask & (((uint64_t)in32[3][33] >>21) + ((uint64_t)in32[3][34] <<8) + ((uint64_t)in32[3][35] <<37));
}
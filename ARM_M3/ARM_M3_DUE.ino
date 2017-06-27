// Name: Julio Santos - julio.ppgca@gmail.com
// Date: nov. 2016
// DMA controling DAC and ADC Sampling.
// Good references to get this working. thanks.
// http://forum.arduino.cc/index.php?topic=224672.0
// https://gist.github.com/pklaus/5921022
// https://forum.arduino.cc/index.php?topic=205096.0
// https://www.keil.com/pack/doc/CMSIS/DSP/html/index.html

#define SAMPLES_NUM   256 // number of point to DAC and ADC converter
#define SAMPLES_AVG   4   // software average of ADC buffer
//#define SINE_WAVE_CLEAN        // comment this to get a noise (6kHz sine) SineWave form
//#define SINE_WAVE_N1 
#define SINE_WAVE_N2

#define ADC_OFFSET    2048 // 12bits adc -> 4096/2=2048
#define ADC_SCALE     (3.3/4095) // max. adc in=3.3V/adc_res-1 -> 3.3V/(4096-1) ~= 0.806mv 

// DAC variables definition
int       freq_inhz  =   60;//1000; // Default Hz
int       freq_intc  =   0; 
volatile uint16_t sptr = 0;

// ADC variables definition
volatile uint16_t bufn,obufn;
uint16_t buf[4][256];   // 4 buffers of 256 readings
    
#ifdef SINE_WAVE_CLEAN 
// t=linspace(0,1/60,2^8);  
// y1 = 0.8*sin(2*pi*60*t);  
// sinal=ceil(y1*2048)+2047;  
// RMS: 0.565 Volts  
uint16_t  SineWave[2][SAMPLES_NUM] = 
{
  {
  2047, 2088, 2128, 2168, 2209, 2249, 2289, 2329, 2368, 2408, 2447, 2486, 2525, 2563, 2602, 2639, 
  2677, 2714, 2751, 2787, 2823, 2858, 2893, 2927, 2961, 2994, 3027, 3059, 3090, 3121, 3151, 3181, 
  3210, 3238, 3265, 3292, 3318, 3343, 3367, 3391, 3413, 3435, 3456, 3476, 3496, 3514, 3532, 3548, 
  3564, 3579, 3593, 3606, 3618, 3629, 3639, 3648, 3656, 3663, 3669, 3675, 3679, 3682, 3684, 3686, 
  3686, 3685, 3683, 3681, 3677, 3672, 3667, 3660, 3652, 3644, 3634, 3623, 3612, 3599, 3586, 3572, 
  3556, 3540, 3523, 3505, 3486, 3466, 3446, 3424, 3402, 3379, 3355, 3330, 3305, 3278, 3251, 3224, 
  3195, 3166, 3136, 3106, 3075, 3043, 3011, 2978, 2944, 2910, 2875, 2840, 2805, 2769, 2732, 2695, 
  2658, 2620, 2583, 2544, 2506, 2467, 2428, 2388, 2349, 2309, 2269, 2229, 2189, 2148, 2108, 2068, 
  }, 
  {
  2027, 1987, 1947, 1906, 1866, 1826, 1786, 1746, 1707, 1667, 1628, 1589, 1551, 1512, 1475, 1437, 
  1400, 1363, 1326, 1290, 1255, 1220, 1185, 1151, 1117, 1084, 1052, 1020, 989,  959,  929,  900,  
  871,  844,  817,  790,  765,  740,  716,  693,  671,  649,  629,  609,  590,  572,  555,  539,  
  523,  509,  496,  483,  472,  461,  451,  443,  435,  428,  423,  418,  414,  412,  410,  409,  
  409,  411,  413,  416,  420,  426,  432,  439,  447,  456,  466,  477,  489,  502,  516,  531,  
  547,  563,  581,  599,  619,  639,  660,  682,  704,  728,  752,  777,  803,  830,  857,  885,  
  914,  944,  974,  1005, 1036, 1068, 1101, 1134, 1168, 1202, 1237, 1272, 1308, 1344, 1381, 1418, 
  1456, 1493, 1532, 1570, 1609, 1648, 1687, 1727, 1766, 1806, 1846, 1886, 1927, 1967, 2007, 2047, 
  }
};
#endif 

#ifdef SINE_WAVE_N1 
// t=linspace(0,1/60,2^8);  
// y2 = 0.5*sin(2*pi*60*t)+0.3*sin(2*pi*120*t)+0.05*sin(2*pi*180*t);  
// sinal=ceil(y2*2048)+2047;  
// RMS: 0.413 Volts  
uint16_t  SineWave[2][SAMPLES_NUM] = 
{
  {
  2047, 2111, 2173, 2236, 2298, 2360, 2421, 2482, 2541, 2600, 2657, 2713, 2768, 2821, 2873, 2923, 
  2971, 3018, 3063, 3105, 3146, 3185, 3221, 3255, 3287, 3317, 3345, 3370, 3392, 3413, 3431, 3447, 
  3460, 3471, 3480, 3486, 3490, 3492, 3492, 3490, 3486, 3479, 3471, 3461, 3449, 3436, 3420, 3404, 
  3385, 3366, 3345, 3322, 3299, 3275, 3249, 3223, 3196, 3168, 3140, 3111, 3082, 3052, 3022, 2992, 
  2962, 2931, 2901, 2871, 2841, 2811, 2781, 2752, 2724, 2695, 2667, 2640, 2613, 2587, 2562, 2537, 
  2512, 2489, 2466, 2444, 2423, 2402, 2382, 2363, 2345, 2328, 2311, 2295, 2279, 2265, 2251, 2238, 
  2225, 2213, 2202, 2191, 2181, 2171, 2162, 2154, 2146, 2139, 2132, 2125, 2119, 2113, 2107, 2102, 
  2097, 2093, 2089, 2085, 2081, 2077, 2074, 2071, 2068, 2065, 2062, 2059, 2056, 2054, 2051, 2049, 
  }, 
  {
  2046, 2044, 2041, 2039, 2036, 2033, 2030, 2027, 2024, 2021, 2018, 2014, 2010, 2006, 2002, 1998, 
  1993, 1988, 1982, 1976, 1970, 1963, 1956, 1949, 1941, 1933, 1924, 1914, 1904, 1893, 1882, 1870, 
  1857, 1844, 1830, 1816, 1800, 1784, 1767, 1750, 1732, 1713, 1693, 1672, 1651, 1629, 1606, 1583, 
  1558, 1533, 1508, 1482, 1455, 1428, 1400, 1371, 1343, 1314, 1284, 1254, 1224, 1194, 1164, 1133, 
  1103, 1073, 1043, 1013, 984,  955,  927,  899,  872,  846,  820,  796,  773,  750,  729,  710,  
  691,  675,  659,  646,  634,  624,  616,  609,  605,  603,  603,  605,  609,  615,  624,  635,  
  648,  664,  682,  703,  725,  750,  778,  808,  840,  874,  910,  949,  990,  1032, 1077, 1124, 
  1172, 1222, 1274, 1327, 1382, 1438, 1495, 1554, 1613, 1674, 1735, 1797, 1859, 1922, 1984, 2047, 
  }
};
#endif

#ifdef SINE_WAVE_N2 
// t=linspace(0,1/60,2^8);  
// y3 = 0.5*sin(2*pi*60*t)+0.1*sin(2*pi*720*t);  
// sinal=ceil(y3*2048)+2047;  
// RMS: 0.360 Volts  
uint16_t  SineWave[2][SAMPLES_NUM] = 
{
  {
  2047, 2132, 2212, 2282, 2338, 2377, 2399, 2403, 2392, 2367, 2335, 2299, 2265, 2238, 2222, 2220, 
  2236, 2269, 2319, 2383, 2458, 2539, 2621, 2699, 2767, 2822, 2862, 2883, 2886, 2872, 2845, 2808, 
  2766, 2724, 2688, 2662, 2649, 2653, 2673, 2711, 2763, 2827, 2898, 2971, 3041, 3102, 3151, 3184, 
  3200, 3197, 3177, 3142, 3096, 3043, 2990, 2940, 2899, 2871, 2858, 2863, 2884, 2920, 2969, 3026, 
  3087, 3145, 3196, 3236, 3261, 3268, 3256, 3227, 3183, 3126, 3062, 2995, 2931, 2874, 2829, 2799, 
  2787, 2791, 2811, 2844, 2887, 2934, 2981, 3022, 3053, 3069, 3069, 3050, 3013, 2960, 2895, 2820, 
  2742, 2666, 2596, 2537, 2492, 2465, 2455, 2461, 2481, 2513, 2550, 2588, 2622, 2647, 2659, 2654, 
  2632, 2592, 2536, 2466, 2386, 2302, 2218, 2140, 2072, 2018, 1981, 1962, 1960, 1973, 1997, 2030, 
  }, 
  {
  2065, 2098, 2122, 2135, 2133, 2114, 2077, 2023, 1955, 1877, 1793, 1709, 1629, 1559, 1503, 1463, 
  1441, 1436, 1448, 1473, 1507, 1545, 1582, 1614, 1634, 1640, 1630, 1603, 1558, 1499, 1429, 1353, 
  1275, 1200, 1135, 1082, 1045, 1026, 1026, 1042, 1073, 1114, 1161, 1208, 1251, 1284, 1304, 1308, 
  1296, 1266, 1221, 1164, 1100, 1033, 969,  912,  868,  839,  827,  834,  859,  899,  950,  1008, 
  1069, 1126, 1175, 1211, 1232, 1237, 1224, 1196, 1155, 1105, 1052, 999,  953,  918,  898,  895,  
  911,  944,  993,  1054, 1124, 1197, 1268, 1332, 1384, 1422, 1442, 1446, 1433, 1407, 1371, 1329, 
  1287, 1250, 1223, 1209, 1212, 1233, 1273, 1328, 1396, 1474, 1556, 1637, 1712, 1776, 1826, 1859, 
  1875, 1873, 1857, 1830, 1796, 1760, 1728, 1703, 1692, 1696, 1718, 1757, 1813, 1883, 1963, 2047, 
  }
};
#endif 


int tcToFreq( int tc_cntr);
int freqToTc( int freq_hz);
void DACC_Handler(void);
void setup_pio_TIOA0();
void TC_setup ();
void DAC_setup ();
void ADC_setup();
void ADC_Handler();

void setup()
{
  Serial.begin (115200); 
  DAC_setup();      
  ADC_setup();  
  freq_intc = freqToTc(freq_inhz);
  TC_setup();        
  setup_pio_TIOA0(); 

  pinMode(12,OUTPUT);
  pinMode(13,OUTPUT);
  digitalWrite(12,LOW);
  digitalWrite(13,LOW);
}

void loop() 
{
  uint16_t g_buff[SAMPLES_NUM/4]={};
  uint16_t sample_id=0;
  uint16_t i,j;
  
    while(true)
    { 
      while(obufn==bufn);    // wait for buffer to be full by DMA
      digitalWrite(12,HIGH); // Osciloscope debug - pin 12 - 3.3V
      for(i=0,j=0;i<SAMPLES_NUM;i+=SAMPLES_AVG) // get 64 aveg. points from 256 sampled buffer
      {
        g_buff[j++] = (buf[obufn][i]   + buf[obufn][i+1] +
                       buf[obufn][i+2] + buf[obufn][i+3] ) / SAMPLES_AVG;    
      }      
      Serial.write(sample_id);
      Serial.write(sample_id>>8);
      Serial.write((uint8_t *)g_buff,128);  // send it - 128 bytes = 64 uint16_t
      sample_id++;
      obufn=(obufn+1)&3; 
      digitalWrite(12,LOW); // Osciloscope debug - pin 12 - 0V (222us)
    }
}

int tcToFreq( int tc_cntr)
{
  if( tc_cntr == 0 ) return 1000;
  return ( 42000000UL / tc_cntr) / (SAMPLES_NUM);   
}

int freqToTc( int freq_hz)
{
  if( freq_hz == 0 ) return 25;
  return ( 42000000UL / freq_hz) / (SAMPLES_NUM);
}

void DACC_Handler(void)
{
  if((dacc_get_interrupt_status(DACC) & DACC_ISR_ENDTX) == DACC_ISR_ENDTX) {
    ++sptr; 
    sptr &=  0x01;

    DACC->DACC_TNPR =  (uint32_t)  SineWave[sptr];      // next DMA buffer
    DACC->DACC_TNCR =  SAMPLES_NUM/2;
  }
}

void setup_pio_TIOA0()  
{
  PIOB->PIO_PDR = PIO_PB25B_TIOA0;  
  PIOB->PIO_IDR = PIO_PB25B_TIOA0;  
  PIOB->PIO_ABSR |= PIO_PB25B_TIOA0;
}


void TC_setup ()
{
  pmc_enable_periph_clk(TC_INTERFACE_ID + 0*3 + 0); 

  TcChannel * t = &(TC0->TC_CHANNEL)[0];            
  t->TC_CCR = TC_CCR_CLKDIS;                        
  t->TC_IDR = 0xFFFFFFFF;                           
  t->TC_SR;                                         
  t->TC_CMR = TC_CMR_TCCLKS_TIMER_CLOCK1 |          
              TC_CMR_WAVE |                         
              TC_CMR_WAVSEL_UP_RC |                 
              TC_CMR_EEVT_XC0 |     
              TC_CMR_ACPA_CLEAR | TC_CMR_ACPC_CLEAR |
              TC_CMR_BCPB_CLEAR | TC_CMR_BCPC_CLEAR;
  
  t->TC_RC = freq_intc;
  t->TC_RA = freq_intc /2;       
  t->TC_CMR = (t->TC_CMR & 0xFFF0FFFF) | TC_CMR_ACPA_CLEAR | TC_CMR_ACPC_SET; 
  t->TC_CCR = TC_CCR_CLKEN | TC_CCR_SWTRG;    
}

void DAC_setup ()
{
  pmc_enable_periph_clk (DACC_INTERFACE_ID) ; // start clocking DAC
  dacc_reset(DACC);
  dacc_set_transfer_mode(DACC, 0);
  dacc_set_power_save(DACC, 0, 1);            // sleep = 0, fastwkup = 1
  dacc_set_analog_control(DACC, DACC_ACR_IBCTLCH0(0x02) | DACC_ACR_IBCTLCH1(0x02) | DACC_ACR_IBCTLDACCORE(0x01));
  dacc_set_trigger(DACC, 1);
  
  dacc_set_channel_selection(DACC, 1);
  dacc_enable_channel(DACC, 1);
  //dacc_set_channel_selection(DACC, 0);
  //dacc_enable_channel(DACC, 0);

  NVIC_DisableIRQ(DACC_IRQn);
  NVIC_ClearPendingIRQ(DACC_IRQn);
  NVIC_EnableIRQ(DACC_IRQn);
  dacc_enable_interrupt(DACC, DACC_IER_ENDTX);

  DACC->DACC_TPR  =  (uint32_t)  SineWave[0];      // DMA buffer
  DACC->DACC_TCR  =  SAMPLES_NUM;
  DACC->DACC_TNPR =  (uint32_t)  SineWave[1];      // next DMA buffer
  DACC->DACC_TNCR =  SAMPLES_NUM/2;
  DACC->DACC_PTCR =  0x00000100;  //TXTEN - 8, RXTEN - 1.
}

void ADC_setup()
{
   pmc_enable_periph_clk(ID_ADC);
  adc_init(ADC, SystemCoreClock, ADC_FREQ_MAX, ADC_STARTUP_FAST);
  ADC->ADC_MR = (ADC->ADC_MR & 0xFFFFFFF0) | (1 << 1) | ADC_MR_TRGEN ;  // 1 = trig source TIO from TC0

  ADC->ADC_CHER=0x80; 

  NVIC_EnableIRQ(ADC_IRQn);
  ADC->ADC_IDR=~(1<<27);
  ADC->ADC_IER=1<<27;
  ADC->ADC_RPR=(uint32_t)buf[0];   // DMA buffer
  ADC->ADC_RCR=256;
  ADC->ADC_RNPR=(uint32_t)buf[1]; // next DMA buffer
  ADC->ADC_RNCR=256;
  bufn=obufn=1;
  ADC->ADC_PTCR=1;
  ADC->ADC_CR=2;
}

/*
    Move DMA pointers to next buffer (transfer completed)
    On Udoo Dual Lite this function takes only 2.9 us to complete.
    It is calles every 16.67us or 60Hz.
*/
void ADC_Handler() 
{    
 digitalWrite(13,HIGH); // Osciloscope debug - pin 13 -> 3.3V
  int f=ADC->ADC_ISR;
  if (f&(1<<27)){
   bufn=(bufn+1)&3;
   ADC->ADC_RNPR=(uint32_t)buf[bufn];
   ADC->ADC_RNCR=256;
  } 
 digitalWrite(13,LOW);  // Osciloscope debug - pin 13 -> 0V
}

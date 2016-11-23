
// Júlio Santos
// gcc serial.c -o serial -lpthread -lm

#include <termios.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>

#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdlib.h>
#include <stdint.h>

#include <math.h>
#include <complex.h>

#define MSG_LENGTH		130 // MSG: IdL,IdH,D[0]L,D[0]H,...,D[63]L,D[64]H
#define BUFFER_SIZE 	MSG_LENGTH*10
#define SERIAL_PORT 	"/dev/ttymxc3"
#define BAUDRATE		B115200
#define PARITY			0	// 0 = none
#define NON_BLOCKING	0

#define SAMPLE_SIZE 	64
#define MIN_POWER		0.06
#define PI 				3.1415926535897932384

#define THREADS_PER_EXEC	1 //10 default
#define WRITE_LOG_FILE
//#define DEBUG_MSG
#define DEBUG_RESULT_VALUES



uint8_t buf[BUFFER_SIZE];
uint16_t g_buff[BUFFER_SIZE]; 		// to do the magic!!! muhahahaha
sem_t s_prod, s_ctr_buf, s_wr_buf;	// buffer ready, buffer is not full, enabled w/r in buffer
char *portname = SERIAL_PORT;
int fd;
time_t mytime=NULL;

void *t_prod(void *none);		// producer function prototype, serial input
void *t_cons(void *none);		// Consumer function prototype
int set_interface_attribs (int fd, int speed, int parity);
void set_blocking (int fd, int should_block);

int log2_r(int N);    /*function to calculate the log2_r(.) of int numbers*/
int check(int n);    //checking if the number of element is a power of 2
int reverse(int N, int n);   //calculating revers number
void ordina(double complex * f1, int N);
void transform(double complex * f, int N);
void FFT(double complex * f, int N, double d);
float RMS(float * f, int n);

int main(int argc, char** argv)
{
	pthread_t tid1, tid2;	// thread identifications
	sem_init(&s_prod,0,0);	// nothing produced yet
	sem_init(&s_ctr_buf,0,BUFFER_SIZE/MSG_LENGTH); //buffer is empty
	sem_init(&s_wr_buf,0,1); //ready to produce
	void *received, *processed;

	int num_exec=1;
	if(argc>1) num_exec=atoi(argv[1]);

	fd = open (portname, O_RDWR | O_NOCTTY | O_SYNC);
	if (fd < 0)
	{
		printf("error %d opening %s: %s", errno, portname, strerror (errno));
		return -1;
	}

	set_interface_attribs (fd, BAUDRATE, PARITY);	// set speed to 115,200 bps, 8n1 (no parity)
	set_blocking (fd, MSG_LENGTH);//NON_BLOCKING);		// set no blocking

	uint32_t i,j;
	for(j=0;j<num_exec;j++)
	{
		for(i=0;i<THREADS_PER_EXEC;i++)
		{
			pthread_create(&tid1, NULL, t_prod, NULL);
			pthread_create(&tid2, NULL, t_cons, NULL);

		}
		for(i=0;i<THREADS_PER_EXEC;i++)
		{
			pthread_join(tid1, &received);
			pthread_join(tid2, &processed);

		}
	}
	printf("\n\n\t\t\tReceived:\t %d messages.",(uint32_t)received);
	printf("\n\t\t\tProcessed:\t %d messages.\n",(uint32_t)processed);

}

void *t_prod(void *none)
{
	int n,i;
	static uint32_t put=0; 	// index to put in the buffer
	static uint32_t received=0;


	sem_wait(&s_ctr_buf); // buffer is not full? Proceed..
	sem_wait(&s_wr_buf); // there is no read operation going on..

	do{
		n = read (fd, buf, sizeof buf);
	}
	while( n < MSG_LENGTH); 	// wait to get all samples

	//printf("\nRecived.: %d bytes.\n",n);
	for (i = 0; i < MSG_LENGTH; i+=2)
	{
		g_buff[put]=buf[i]|(buf[i+1]<<8);
		if(++put>=BUFFER_SIZE) put=0; // circular type
	}
	received++;
#ifdef DEBUG_MSG
	//printf("\n\t\t\tRecived.: %d bytes.\n",n);
	printf("Get sample ID: %d\n",buf[0]|buf[1]<<8);
#endif

	tcflush(fd,TCIOFLUSH);
	sem_post(&s_wr_buf); // release buffer w/r operation
	sem_post(&s_prod); 	// item produced and buffered

	return (void *)received;//NULL;
}

void *t_cons(void *none)
{
	static uint32_t get=0; // index to get from the buffer
	uint32_t i;
	uint16_t cons=0;
	static uint32_t processed;

	uint16_t signal[SAMPLE_SIZE];
	uint16_t sample_id;
	float rms=0,power=0;
	double complex vec[SAMPLE_SIZE];

	//usleep(500000);

	sem_wait(&s_prod); // is there an item produced? Proceed.
	sem_wait(&s_wr_buf); // block buffer w/r

	sample_id=g_buff[get++];
	for(i=0;i<SAMPLE_SIZE;i++)
	{
		signal[i]=g_buff[get];
		if(++get>=BUFFER_SIZE) get=0;
	}
	processed++;
#ifdef DEBUG_MSG
	printf("\t\t\t\t\t\t Set sample ID: %d\n",sample_id);
#endif


#ifdef DEBUG_RESULT_VALUES

	// RMS
	for(i=0;i<SAMPLE_SIZE;i++)
		rms+=(((float)signal[i]-2048)*(3.3/4095))*(((float)signal[i]-2048)*(3.3/4095));
	rms=sqrt(rms/SAMPLE_SIZE);
	printf("\nRMS[%d]: %.3f \n",sample_id,rms);

	//	FFT
	for(i=0;i<SAMPLE_SIZE;i++) vec[i]=((float)signal[i]-2048)*(3.3/4095);
	FFT(vec, SAMPLE_SIZE, 1);
	printf("Looking for frequencies greater than %.2f dB\n",MIN_POWER);
	for(i = 0; i < SAMPLE_SIZE/2; i++)
	{
		power=2*cabs(vec[i])/SAMPLE_SIZE;
		if(power>MIN_POWER) printf("Freq. %dHz \tPower: %.4lf dB\n",i*60,power);
	}

#endif
	// data file
	FILE *f1 = fopen("data.csv", "w");
	if (f1 == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	for(i=0;i<SAMPLE_SIZE-1;i++) fprintf(f1, "%d\n", signal[i]);
	fprintf(f1,"%d",signal[i++]);
	fclose(f1);
#ifdef WRITE_LOG_FILE

	// log file
	FILE *f2 = fopen("log.txt", "a");
	if (f2 == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	time(&mytime);
	fprintf(f2, "\n%s", ctime(&mytime));
	fprintf(f2, "Sample number:\t%d\n", sample_id);
	fprintf(f2,"RMS Value: %.3f Volts\n",sqrtf(rms));
	fprintf(f2,"FFT Analyses:\nLooking for frequencies greater than %.2f dB\n",MIN_POWER);
	for(i=0;i<SAMPLE_SIZE/2;i++)
	{
		power=2*cabs(vec[i])/SAMPLE_SIZE;
		if(power>MIN_POWER) fprintf(f2,"\tFreq. %dHz \tPower: %.4lf dB\n",i*60,power);
	}
	fclose(f2);
#endif

	sem_post(&s_wr_buf); // realease buffer w/r
	sem_post(&s_ctr_buf); // item consumed, free buffer position

	return (void *)processed;
}

int set_interface_attribs (int fd, int speed, int parity)
{
	struct termios tty;
	memset (&tty, 0, sizeof tty);
	if (tcgetattr (fd, &tty) != 0)
	{
		printf ("error %d from tcgetattr", errno);
		return -1;
	}

	cfsetospeed (&tty, speed);
	cfsetispeed (&tty, speed);

	tty.c_cflag = (tty.c_cflag & ~CSIZE) | CS8;     // 8-bit chars
	// disable IGNBRK for mismatched speed tests; otherwise receive break
	// as \000 chars
	tty.c_iflag &= ~IGNBRK;         // ignore break signal
	tty.c_lflag = 0;                // no signaling chars, no echo,
	// no canonical processing
	tty.c_oflag = 0;                // no remapping, no delays
	tty.c_cc[VMIN]  = 0;            // read doesn't block
	tty.c_cc[VTIME] = 5;            // 0.5 seconds read timeout

	tty.c_iflag &= ~(IXON | IXOFF | IXANY); // shut off xon/xoff ctrl

	tty.c_cflag |= (CLOCAL | CREAD);// ignore modem controls,
	// enable reading
	tty.c_cflag &= ~(PARENB | PARODD);      // shut off parity
	tty.c_cflag |= parity;
	tty.c_cflag &= ~CSTOPB;
	tty.c_cflag &= ~CRTSCTS;

	if (tcsetattr (fd, TCSANOW, &tty) != 0)
	{
		printf ("error %d from tcsetattr", errno);
		return -1;
	}
	return 0;
}

void set_blocking (int fd, int should_block)
{
	struct termios tty;
	memset (&tty, 0, sizeof tty);
	if (tcgetattr (fd, &tty) != 0)
	{
		printf ("error %d from tggetattr", errno);
		return;
	}

	tty.c_cc[VMIN]  = should_block ? 1 : 0;
	tty.c_cc[VTIME] = 5;            	// 0.5 seconds read timeout

	tcflush(fd, TCIFLUSH);
	if (tcsetattr (fd, TCSANOW, &tty) != 0)
		printf ("error %d setting term attributes", errno);
}

/*** Fourier ***/
int log2_r(int N)    /*function to calculate the log2_r(.) of int numbers*/
{
	int k = N, i = 0;
	while(k) {
		k >>= 1;
		i++;
	}
	return i - 1;
}

int check(int n)    //checking if the number of element is a power of 2
{
	return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n)    //calculating revers number
{
	int j, p = 0;
	for(j = 1; j <= log2_r(N); j++) {
		if(n & (1 << (log2_r(N) - j)))
			p |= 1 << (j - 1);
	}
	return p;
}

void ordina(double complex * f1, int N) //using the reverse order in the array
{
	int i,j;
	double complex f2[SAMPLE_SIZE];
	for(i = 0; i < N; i++)
		f2[i] = f1[reverse(N, i)];
	for(j = 0; j < N; j++)
		f1[j] = f2[j];
}

void transform(double complex * f, int N) //
{
	int i,j;
	ordina(f, N);           //first: reverse order
	double complex * W;
	W = (double complex *)malloc(N / 2 * sizeof(double complex));
	W[1] = 1 + (-2. * PI / N)*I;
	W[0] = 1;
	for(i = 2; i < N / 2; i++)
		W[i] = cpow(W[1], i);
	int n = 1;
	int a = N / 2;
	for(j = 0; j < log2_r(N); j++) {
		for(i = 0; i < N; i++) {
			if(!(i & n)) {
				double complex temp = f[i];
				double complex Temp = W[(i * a) % (n * a)] * f[i + n];
				f[i] = temp + Temp;
				f[i + n] = temp - Temp;
			}
		}
		n *= 2;
		a /= 2;
	}
}

void FFT(double complex * f, int N, double d)
{
	int i;
	transform(f, N);
	for(i = 0; i < N; i++)
		f[i] *= d; //multiplying by step
}

float RMS(float * f, int n)
{
	int i;
	float rms=0;
	for(i=0;i<n;i++) rms+=(f[i]*f[i]);
	return sqrtf(rms/n);
}
/*** end Fourier ***/



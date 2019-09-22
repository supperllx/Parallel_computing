#include <stdio.h>
#include <pthread.h>
#include <iostream>

#include "CycleTimer.h"

using namespace std;

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);


//
// workerThreadStart --
//
// Thread entrypoint.
void* workerThreadStart(void* threadArgs) {

    WorkerArgs* args = static_cast<WorkerArgs*>(threadArgs);

    // TODO: Implement worker thread here.
    double minSerial = 1e30;
    double startTime,endTime;
    startTime = CycleTimer::currentSeconds();

    float x0 = args->x0;
    float y0 = args->y0;
    float x1 = args->x1;
    float y1 = args->y1;
    int width = args->width;
    int height = args->height;
    int maxIterations = args->maxIterations;
    
    int num_rows = 1;
    int start_row;
    for(int i=0;i!=height/args->numThreads;++i){
        start_row = args->threadId + i*args->numThreads;
        mandelbrotSerial(x0,y0,x1,y1,width, height, start_row, num_rows,maxIterations,args->output);
    }
    /*
    int num_rows = height/args->numThreads;
    int start_row = num_rows * args->threadId;

    mandelbrotSerial(x0, y0, x1, y1, width, height, start_row, num_rows, maxIterations, args->output);
    */
    endTime = CycleTimer::currentSeconds();
    minSerial = std::min(minSerial, endTime - startTime);
    cout<<"The thread: "<<args->threadId<<" use: "<<minSerial<<"s"<<endl;

    return NULL;
}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Multi-threading performed via pthreads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    const static int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    pthread_t workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i=0; i<numThreads; i++) {
        // TODO: Set thread arguments here.
        args[i].numThreads = numThreads;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;

        args[i].x0 = x0;
        //args[i].y0 = i*(height/numThreads)+y0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        //args[i].y1 = i*(height/numThreads)+y1;
        args[i].y1 = y1;

        args[i].output = output;
        args[i].threadId = i;
    }

    // Fire up the worker threads.  Note that numThreads-1 pthreads
    // are created and the main app thread is used as a worker as
    // well.

    int ret;

    for (int i=0; i<numThreads; i++){
        ret = pthread_create(&workers[i], NULL, workerThreadStart, &args[i]);
        if(ret){
            cout<<"error in creating thread: "<<i<<endl;
            exit(-1);
        }
    }

    // wait for worker threads to complete
    for (int i=0; i<numThreads; i++)
        pthread_join(workers[i], NULL);
}

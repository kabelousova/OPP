#include <iostream>
#include <mpi.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <climits>
#include <cstddef>
#include <vector>

#define COUNT_TASKS 10000
#define COUNT_ITER 1000
#define DEFAULT_WEIGHT 100
#define EXTRA_WEIGHT 500
#define REQUEST_TAG 101
#define ANSWER_TAG 102
#define DATA_TAG 103
#define AN_T 10

struct TaskInfo {
    double start;
    double step;
    int repeatNum;
};

pthread_mutex_t taskMutex;
std::vector <TaskInfo> tasks;
MPI_Datatype TASK_INFO_TYPE;

int currentRank;
int countRanks;
int iterCounter = 0;
int currentTask;
int answerTasks;
double result = 0.f;

int getNewTasks(int procRank) {
    int req = 1;

    MPI_Send(&req, 1, MPI_INT, procRank, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&req, 1, MPI_INT, procRank, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (req == 0) {
        return 1;
    }

    TaskInfo task[answerTasks];
    MPI_Recv(&task, answerTasks, TASK_INFO_TYPE, procRank, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    pthread_mutex_lock(&taskMutex);

    for (size_t i = 0; i < answerTasks; i++) {
        tasks.push_back(task[i]);
    }

    pthread_mutex_unlock(&taskMutex);
    return 0;
}

double randForTask(double a, double b) {
    double result = (double) rand() / INT_MAX;
    return result * (b - a) + a;
}

void generateTasks(int countTasks) {
    int weight = abs(currentRank - (iterCounter % countRanks));
    tasks.clear();

    for (int i = 0; i < countTasks; ++i) {
        TaskInfo task;
        task.start = randForTask(10000., 50000.);
        task.repeatNum = DEFAULT_WEIGHT + EXTRA_WEIGHT * weight;
        task.step = 1 / randForTask(2., 1000.);
        if (task.step >= 1) {
            task.step = 1 / tasks[i].step;
        }
        tasks.push_back(task);
    }
}


double doTask(TaskInfo &task) {
    double sum = 0;
    double tmp = task.start;

    for (int i = 0; i < task.repeatNum; ++i) {
        tmp = sqrt(tmp);
        sum += tmp;
        tmp *= task.
        step;
    }

    return sum;
}

void mainThread() {
    int *processReqested = new
    int[countRanks];
    double start_time = MPI_Wtime();
    for (iterCounter = 0; iterCounter < COUNT_ITER; ++iterCounter) {
        pthread_mutex_lock(&taskMutex);
        int count = COUNT_TASKS / countRanks;
        if (currentRank <= countRanks - (COUNT_TASKS % countRanks)) {
            ++count;
        }
        generateTasks(count);
        pthread_mutex_unlock(&taskMutex);

        for (int i = 0; i < countRanks; ++i) {
            processReqested[i] = 0;
        }
        processReqested[currentRank] = 1;

        pthread_mutex_lock(&taskMutex);
        currentTask = 0;
        for (int i = 0; i < countRanks;) {
            while (currentTask < tasks.size()) {
                ++currentTask;
                pthread_mutex_unlock(&taskMutex);
                result += doTask(tasks[currentTask - 1]);
                pthread_mutex_lock(&taskMutex);
            }
            pthread_mutex_unlock(&taskMutex);

            int procRank = (currentRank + i) % countRanks;
            if (processReqested[procRank]) {
                ++i;
                continue;
            }
            processReqested[procRank] = getNewTasks(procRank);
            if (processReqested[procRank]) {
                ++i;
            }

            pthread_mutex_lock(&taskMutex);
        }
        double end_time = MPI_Wtime();
        printf("Current iter: %d, tasks complete: %d, global res: %lf, time spent: %lf", iterCounter, currentTask,
               result,
               end_time - start_time);
        start_time = end_time;
        pthread_mutex_unlock(&taskMutex);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int req = 0;
    MPI_Send(&req, 1, MPI_INT, currentRank, REQUEST_TAG, MPI_COMM_WORLD);
    delete[]
    processReqested;
}

void *dataThread(void *args) {
    while (iterCounter < COUNT_ITER) {
        MPI_Status status;
        int answer;
        int res;

        MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);

        if (!res) {
            break;
        }

        pthread_mutex_lock(&taskMutex);
        if (currentTask >= tasks.size() - answerTasks) {
            answer = 0;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            pthread_mutex_unlock(&taskMutex);
            continue;
        }

        answer = 1;
        MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
        MPI_Send(&tasks.back() - answerTasks - 1, answerTasks, TASK_INFO_TYPE, status.MPI_SOURCE, DATA_TAG,
                 MPI_COMM_WORLD);

        for (size_t i = 0; i < answerTasks; i++) {
            tasks.pop_back();
        }

        pthread_mutex_unlock(&taskMutex);
    }

    return nullptr;
}

void generate() {
    pthread_t thread;
    pthread_create(&thread, nullptr, dataThread, nullptr);

    mainThread();

    pthread_join(thread, nullptr);
}

void createType() {
    int blocklengths[] = {1, 1, 1};

    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(TaskInfo, start);
    offsets[1] = offsetof(TaskInfo, step);
    offsets[2] = offsetof(TaskInfo, repeatNum);

    MPI_Type_create_struct(3, blocklengths, offsets, types, &TASK_INFO_TYPE);
    MPI_Type_commit(&TASK_INFO_TYPE);
}

int main(int argc, char **argv) {
    double start, end;
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &countRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);

    srand(time(nullptr));

    answerTasks = AN_T;
    if (!currentRank)
        printf("answerTasks: %d\n", answerTasks);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Can't get need level\n");
        MPI_Finalize();
        return 1;
    }

    createType();

    pthread_mutex_init(&taskMutex, nullptr);

    start = MPI_Wtime();
    generate();
    end = MPI_Wtime();

    if (!currentRank)
        printf("Time: %lf\n", end - start);

    pthread_mutex_destroy(&taskMutex);

    MPI_Finalize();

    return 0;
}

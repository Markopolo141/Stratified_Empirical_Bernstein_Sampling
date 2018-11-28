

void* salloc(const char *pathname, int proj_id, size_t size, int shmflg) {
	key_t key = ftok(pathname, proj_id);
	int shmid = shmget(key, size, shmflg);
	if (shmid == -1) {
		perror("shmget");
		exit(1);
	}
	void *data = shmat(shmid, (void *)0, 0);
	if (data == (char *)(-1)) {
		perror("shmat");
		exit(1);
	}
	if (shmctl(shmid, IPC_RMID, NULL) == -1) {
		perror("shmctl");
		exit(1);
	}
	return data;
}
void sfree(void *shmaddr) {
	if (shmdt(shmaddr) == -1) {
		perror("shmdt");
		exit(1);
	}
}

void delay(long iters) {
	while (iters > 0) {
		iters--;
	}
}


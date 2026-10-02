#!/Local/bfe_maruvka/anaconda3/bin/python -u
import subprocess, psutil, time, sys


def parse_ram_request(ram_request: str) -> int:
    ram_request = ram_request.rstrip().upper()
    if ram_request[-2:]=="MB":
        return int(ram_request[:-2])*10**6
    elif ram_request[-2:]=="GB":
        return int(ram_request[:-2])*10**9
    else:
        raise RuntimeError("Improper memory request format")


class MemoryBlock:
    def __init__(self):
        self.block = None

    @property
    def size(self):
        if self.block is None:
            return 0
        else:
            return len(self.block)

    def resize(self, new_size):
        if new_size == 0:
            self.block = None # gc will get memory back
        else:
            self.block = bytearray(new_size)


# Start subprocess
def main():
    if len(sys.argv) < 3 or sys.argv[1] == "--help":
        print("USAGE: python memory_hog.py N[MB|GB] executable arg1 arg2...\n ex. python memory_hog.py 5GB samtools sort mybam.bam")
        exit()
    desired_ram_str = sys.argv[1]
    desired_ram = parse_ram_request(desired_ram_str)
    # print(sys.argv[2:])
    p = subprocess.Popen(sys.argv[2:])

    proc = psutil.Process(p.pid)
    memory_block = MemoryBlock()

    try:
        while p.poll() is None:  # while still running
            mem = proc.memory_info().rss
            if (mem+memory_block.size) < desired_ram:
                memory_block.resize(desired_ram-(mem+memory_block.size))
            elif (mem+memory_block.size) > int(desired_ram*2):
                memory_block.resize(0)

            # mem = proc.memory_info().rss / (1024 * 1024)  # in MB
            # print(f"Memory usage: {mem:.2f} MB")
            time.sleep(1)
    except psutil.NoSuchProcess:
        raise RuntimeError("ERROR: LOST PROCESSS PID???")

    # print("Process finished.")


if __name__ == '__main__':
    # print(sys.argv)
    # subprocess.Popen(['Start-Sleep', '-Seconds', '10'])
    # subprocess.Popen(['sleep', '10'])
    main()

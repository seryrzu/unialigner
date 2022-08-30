import argparse
import matplotlib
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    params = parser.parse_args()
    min_rare_len = []
    with open(params.input) as f:
        f.readline()
        for line in f:
            min_rare_len.append(int(line.strip().split()[1]))
    min_rare_len = min_rare_len[:len(min_rare_len) // 2]
    plt.plot(min_rare_len)
    plt.ylim(0, 30000)
    plt.title("Length of unique substring starting at each position")
    plt.xlabel("position")
    plt.ylabel("length")
    plt.savefig(params.output, format='pdf')
    plt.close()

if __name__ == "__main__":
    main()
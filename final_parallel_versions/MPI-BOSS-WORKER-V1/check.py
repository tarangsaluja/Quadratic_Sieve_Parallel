"""
This program checks the correctness of our QS implementation, based
on the output of step1.cpp and step2.cpp

Tai Thongthai and Tarang Saluja
"""
def main():

    smooth_nums = []
    factorbase = []
    expo_matrix = []

    with open("Smooth_Num.txt") as sm_file:
        while True:
            line = sm_file.readline()
            if line == "":
                break
            smooth_nums.append(int(line))
    sm_file.close()

    with open("factorbase.txt") as fb_file:
        while True:
            line = fb_file.readline()
            if line == "":
                break
            factorbase.append(int(line.split()[0]))
    fb_file.close


    with open("Expo_Matrix.txt") as em_file:
        while True:
            temp = []
            line = em_file.readline()
            if line == "":
                break
            for power in line:
                if power != "\n":
                    temp.append(int(power))
            expo_matrix.append(temp)

    good = True
    errors = []
    for i in range(len(smooth_nums)):
        res = check(smooth_nums[i], expo_matrix[i], factorbase)
        if not res:
            errors.append("Error at"+ str(smooth_nums[i]))
            good = False

    if good == True:
        print("All good!")
    else:
        print("Errors at: ")
        print(errors)

"""
This function checks whether a row in our exponential matrix
correctly calculates its corresponding smooth number, given the factorbase

Input: Takes in a single smooth number, a relation (exponential matrix row),
and our factorbase 

Returns true if it is correct. False otherwise
"""
def check(smooth, relation, factorbase):

    result = 1
    for i in range(len(relation)):
        result *= factorbase[i]**relation[i]

    if result == smooth:
        return True
    else:
        return False

if __name__ == "__main__":
    main()

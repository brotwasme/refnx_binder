#make fake data
step_in_x = 0.05
no_data_points = 20
sld = 0.2
error = 0.01

def fake_data():
    x,y,z=[],[],[]
    current_x = 0
    for i in range(no_data_points):
        x.append(current_x)
        y.append(sld)
        z.append(error)
        current_x += step_in_x
    return x,y,z

if __name__ == "__main__": #useful for checks in terminal
    print(fake_data())
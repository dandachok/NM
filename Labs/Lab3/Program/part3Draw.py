import matplotlib.pyplot as plt

x = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
y = [-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.313813]
y1 = [-1.312, -0.1901, 0.9315, 2.053, 3.175, 4.296]
y2 = [0.06011, -0.4645, -0.166, 0.9556, 2.9, 5.668]

def draw(x, y, y1, y2):
    fig = plt.figure()
    plt.scatter(x, y)
    plt.plot(x, y1, label = u'First')
    plt.plot(x, y2, label = u'Second')
    plt.grid(True)
    plt.legend()
    plt.show()

#set terminal gif animate delay 100
#set output 'result.gif'

time_step=5e-5

do for [i=0:50] {
    set title sprintf('time=%fs', i*time_step)
    # 2D
    plot sprintf('res_%d.dat', i) w image notitle
    # 3D
    #sp sprintf('res_%d.dat', i) u 1:2:3:4 w p pt 7 ps 0.4 linecolor palette notitle
    pause 0.05
}
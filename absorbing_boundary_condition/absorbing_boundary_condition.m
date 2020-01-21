%absorbing_boundary_condition
function absorbing_boundary_condition %main function
l=1.5;%distance
m=15;%num_dist_interval
T=4;%time
N=40;%num_time_interval
alpha=1;%velocity
h=l/m;%dist_interval,1.5/15=0.1
k=T/N;%time_interval,4/40=0.1
lambda=k*alpha/h;%0.1*1/0.1=1
W=matrix(m,N);
W=step4(W,m,h,k,lambda,@f2,@g);%@call function, f1 for non-noise f2 for noise
W=step5(W,m,N,lambda,h,alpha,k);
draw(W,m,N,T,l);
output(m,N,W);%output table
end

function W=matrix(m,N)%set grid, column=dist, row=time
W=zeros(m+1,N+1);

for i=1:1:N+1
    W(1,i)=0;%dist=0 and 1.5 are zero
    W(m+1,i)=0;
end

W(1,2)=2;%initail condition
W(1,3)=10;
W(1,4)=8;
W(1,5)=5;

end
function x1=step4(W,m,h,k,lambda,f1,g)%setW(j,1) and W(j,2)
for j=2:1:m
   W(j,1)=f1((j-1)*h);%j=grid number, (j-1)*h=grid value
   W(j,2)=(1-lambda^2)*f1((j-1)*h)+(lambda^2)/2*(f1((j-1+1)*h)+f1((j-1-1)*h))+k*g((j-1)*h);%g=partial u/partial t¡Ag=(x,0)=0¡Awhen 0<=x<=1.5
end
x1=W;
end
function x2=step5(W,m,N,lambda,h,alpha,k)%calculate remian gird¡Aexcept for W(1,i)¡AW(16,i)¡AW(j,1)¡AW(j,2)
for i=2:1:N
    for j=1:1:m+1
        if j==1 && i>6
           W(j,i+1)=(1/h)*(alpha*k*W(j+1,i)-(alpha*k-h)*W(j,i));%left-absorbing_boundary_condition
        elseif j==m+1
           W(j,i+1)=(1/h)*(alpha*k*W(j-1,i)-(alpha*k-h)*W(j,i));
        elseif j==1   && i<=6
            %skip
        else
           W(j,i+1)=2*(1-lambda^2)*W(j,i)+(lambda^2)*(W(j+1,i)+W(j-1,i))-W(j,i-1);%right-absorbing_boundary_condition
        end
     end
end
x2=W;
%step 6 
end
function draw(W,m,N,T,l)%plot
x=zeros(m+1,1);%x axis
for j=1:1:m+1
     x(j,1)=0.1*(j-1);
end

t_num=N+1;
x_num=m+1;
u_min=min(min(W));
u_max=max(max(W));
%plot(W);

for t=1:1:t_num
    time=((t_num-t)*0+(t-1)*T)/(t_num-1);%t1=0, t2=T
    plot(x(1:x_num),W(1:x_num,t),'b*-')%plot¡A'b*-' b for blue,  * for point shape,  - for line style
    grid on%show grid
    axis([0,l,u_min-1.0,u_max+1.0]);%x_start, x_end, y_start, y_end
    title(sprintf('step %d,Time %f \n',t-1,time))%title
    xlabel('<---x--->')
    ylabel('vertical displacement')
    pause
end

%make gif
fig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'absorbing_boundary_condition2.gif';
for t=1:1:t_num
    % Draw plot
    time=((t_num-t)*0+(t-1)*T)/(t_num-1);%t1=0, t2=T
    plot(x(1:x_num),W(1:x_num,t),'b*-')%plot¡A'b*-' b for blue,  * for point shape,  - for line style
    grid on%show grid
    axis([0,l,u_min-1.0,u_max+1.0]);%x_start, x_end, y_start, y_end
    title(sprintf('step %d,Time %f \n',t-1,time))%title
    xlabel('<---x--->')
    ylabel('vertical displacement')
    drawnow 
    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.3); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',0.3); 
    end 
end
end
function f1=f1(x)%modeling without background waves
f1=0*sin(pi*x);
end
function f2=f2(x)%modeling with background waves
f2=sin(pi*x);
end
function g1=g(x)
g1=0*(x-1);%g=partial u/partial t¡Ag=(x,0)=0¡Awhen 0<=x<=1.5
end
function output(m,N,W)
fid=fopen('output.txt','W');
fprintf(fid,'output W \r\n');
for j=1:1:m+1
    for i=1:1:N+1
        fprintf(fid,'%d ',W(j,i));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);
end
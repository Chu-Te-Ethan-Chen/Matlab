%shw1 
%f(x)=u(x,0)=sin(pi*x) 0<=x<=1
%step 1 set initial condition
l=1;%distance
m=10;%num_dist_interval
T=1;%time
N=20;%num_time_interval
alpha=2;

h=l/m;%dist_interval,1/10=0.1
k=T/N;%time_interval,1/20=0.05
lambda=k*alpha/h;%0.05*2/0.1=1

%step 2 %grid and initial condition
W=zeros(11,21);
for i=1:1:11
    W(1,i)=0;%first and last point are zero
    W(11,i)=0;
end

%step 3 initial condition
W(1,1)=sin(pi*0);%W(1,1)=f(0)
W(11,1)=sin(pi*1);%W(11,1)=f(l)

%step 4
for j=2:1:10
   W(j,1)=sin(pi*(j-1)*h);
   W(j,2)=(1-lambda^2)*sin(pi*(j-1)*h)+(lambda^2)/2*(sin(pi*(j-1+1)*h)+sin(pi*(j-1-1)*h))+k*0;%g(x)=f(x) partial differential
end

%step 5
for i=2:1:20
    for j=2:1:10
       W(j,i+1)=2*(1-lambda^2)*W(j,i)+(lambda^2)*(W(j+1,i)+W(j-1,i))-W(j,i-1);
    end
end

%step 6 

%output grid
theory=zeros(11,1);
for i=1:1:11
   theory(i,1)=sin(pi*(i-1)*0.1)*cos(2*pi*14*k);
end

fid=fopen('shw1.txt','w');
fprintf(fid,'x            W(x,14)       theory       residual \r\n');
for i=1:1:11
   fprintf(fid,'%f     %f     %f     %e \r\n',(i-1)*0.1,W(i,15),theory(i,1),W(i,15)-theory(i,1));
end
fclose(fid);

%plot
x=zeros(11,1);%x axis
for j=1:1:11
x(j,1)=0.1*(j-1);
end

t_num=N+1;
x_num=m+1;
u_min=min(min(W));
u_max=max(max(W));
plot(W);

for t=1:1:t_num
    time=((t_num-t)*0+(t-1)*1)/(t_num-1);%t1=0
    plot(x(1:x_num),W(1:x_num,t),'b*-')
    grid on
    axis([0,1,u_min-1.0,u_max+1.0]);
    title(sprintf('step %d,Time %f \n',t-1,time))
    xlabel('<--x-->')
    ylabel('vertical displacement')
    %pause    
end 

%make gif
fig = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Amplitude.gif';
for t=1:1:t_num
    % Draw plot
    time=((t_num-t)*0+(t-1)*1)/(t_num-1);
    plot(x(1:x_num),W(1:x_num,t),'b*-')
    grid on
    axis([0,1,u_min-1.0,u_max+1.0]);
    title(sprintf('step %d,Time %f \n',t-1,time))
    xlabel('<--x-->')
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
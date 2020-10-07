%Charged particle dynamics

m_p = 100;            %mass of particles in kg
Np = 4;             %number particles
k = 10;           %proportionality constant Newton's law
L = 50;             %size domain
writevideo = true;  %write to video guard
max_t = 6;          %max simulation time
dt = 0.001;         %time step size

%Initialization
particle_track_1_x = zeros(10*Np,1);
particle_track_1_y = zeros(10*Np,1);
particle_track_2_x = zeros(10*Np,1);
particle_track_2_y = zeros(10*Np,1);
particle_track_3_x = zeros(10*Np,1);
particle_track_3_y = zeros(10*Np,1);

err = 1e-8;
v_x_p = zeros(Np,1);
v_y_p = zeros(Np,1);

v_x_p(2) = 0;
v_y_p(2) = 1200;

v_x_p(3) = 0;
v_y_p(3) = 800;

v_x_p(4) = 0;
v_y_p(4) = 625;

%v_x_p(5) = 300;
%v_y_p(5) = 0;

m_p_i = zeros(Np,1);
for i = 1:Np
    m_p_i(i) = m_p;
end

m_p_i(1) = 1000000;

min_y = -L;
max_y = L;
min_x = -L;
max_x = L;

East_wall = 1;
West_wall = 2;
North_wall = 3;
South_wall = 4;


%Initalize positions of particles
X_p = 2*L*rand(Np,1) - L;
Y_p = 2*L*rand(Np,1) - L;

X_p(1) = 0;
Y_p(1) = 0;

X_p(2) = 10;
Y_p(2) = 0;

X_p(3) = 20;
Y_p(3) = 0;

X_p(4) = 30;
Y_p(4) = 0;

%X_p(5) = 0;
%Y_p(5) = -10;

%Generate video of simulation
frameps = 48;
if writevideo==true
    writerObj = VideoWriter('C:\Users\d-w-h\Desktop\Home\Particle_dynamics.avi','Motion JPEG AVI');
    writerObj.FrameRate = frameps;
    open(writerObj);
end

%Start simulation loop
t = 0;
frame_counter = 0;
timestep = 1;
while t < max_t
    %Check collision times with boundary
    coll_time = 1e+30;
    coll_partner_1 = 0;
    coll_partner_2 = 0;
    collision_with_p = false;
    
    %Checking collision time between walls and positively charged particle
    for n = 1:Np   
        coll_time_east = (L - X_p(n)) / v_x_p(n);
        if coll_time >= coll_time_east && coll_time_east > 0
            coll_time = coll_time_east;
            coll_partner_1 = n;
            coll_partner_2 = East_wall;
            collision_with_p = true;
        end
        
        coll_time_west = (-L - X_p(n)) / v_x_p(n);
        if coll_time >= coll_time_west && coll_time_west > 0
            coll_time = coll_time_west;
            coll_partner_1 = n;
            coll_partner_2 = West_wall;
            collision_with_p = true;
        end 
        
        coll_time_north = (L - Y_p(n)) / v_y_p(n);
        if coll_time >= coll_time_north && coll_time_north > 0
            coll_time = coll_time_north;
            coll_partner_1 = n;
            coll_partner_2 = North_wall;
            collision_with_p = true;            
        end  

        coll_time_south = (-L - Y_p(n)) / v_y_p(n);
        if coll_time >= coll_time_south && coll_time_south > 0
            coll_time = coll_time_south;
            coll_partner_1 = n;
            coll_partner_2 = South_wall;
            collision_with_p = true;            
        end            
    end
    
    %Updating particle velocities and positions
    if coll_time < dt
        X_p = v_x_p * coll_time * (1 - err) + X_p;
        Y_p = v_y_p * coll_time * (1 - err) + Y_p;
                
        t = t + coll_time
                
        %Update velocities after collision with wall
        if coll_partner_2 == East_wall || coll_partner_2 == West_wall
            v_x_p(coll_partner_1) = -v_x_p(coll_partner_1);
        elseif coll_partner_2 == North_wall || coll_partner_2 == South_wall
            v_y_p(coll_partner_1) = -v_y_p(coll_partner_1);
        end

        %Calculate force acting on particles
        Fp = zeros(2, Np);        
        for n_p = 1:Np
            %Calculating force on particles due to other particles
            for n_p_n = 1:Np
                if n_p_n ~= n_p
                    rp = [X_p(n_p); Y_p(n_p)];
                    rp_n = [X_p(n_p_n); Y_p(n_p_n)];
                    r = -(rp - rp_n);
                    Fp(:,n_p) = Fp(:,n_p) + k * m_p_i(n_p) * m_p_i(n_p_n) * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end

        %Updating velocities of p particles
        for n_p = 1:Np
            v_x_p(n_p) = Fp(1,n_p) / m_p_i(n_p) * coll_time + v_x_p(n_p);
            v_y_p(n_p) = Fp(2,n_p) / m_p_i(n_p) * coll_time + v_y_p(n_p);
        end        
        
    elseif coll_time >= dt
        X_p = v_x_p * dt + X_p;
        Y_p = v_y_p * dt + Y_p;
        
        t = t + dt
        
        %Calculate force acting particles
        Fp = zeros(2, Np);        
        for n_p = 1:Np
            %Calculating force on particles due to other particles
            for n_p_n = 1:Np
                if n_p_n ~= n_p
                    rp = [X_p(n_p); Y_p(n_p)];
                    rp_n = [X_p(n_p_n); Y_p(n_p_n)];
                    r = -(rp - rp_n);
                    Fp(:,n_p) = Fp(:,n_p) + k * m_p_i(n_p) * m_p_i(n_p_n) * r / (r'*r * sqrt(r'*r) + err);
                end
            end            
        end
        
        %Updating velocities of particles
        for n_p = 1:Np
            v_x_p(n_p) = Fp(1,n_p) / m_p_i(n_p) * dt + v_x_p(n_p);
            v_y_p(n_p) = Fp(2,n_p) / m_p_i(n_p) * dt + v_y_p(n_p);
        end         
    end
    
    particle_track_1_x(timestep) = X_p(2);
    particle_track_1_y(timestep) = Y_p(2);
    particle_track_2_x(timestep) = X_p(3);
    particle_track_2_y(timestep) = Y_p(3);
    particle_track_3_x(timestep) = X_p(4);
    particle_track_3_y(timestep) = Y_p(4);
    
    %Capture frames for video creation
    if frame_counter == floor(t/(2*dt))
        frame_counter = frame_counter + 1;  
        
        plot(X_p(1), Y_p(1), 'r.', 'MarkerSize', 50)
        hold on
        plot(particle_track_1_x(1:timestep), particle_track_1_y(1:timestep), 'k');
        %plot(particle_track_2_x(1:timestep), particle_track_2_y(1:timestep), 'k');
        %plot(particle_track_3_x(1:timestep), particle_track_3_y(1:timestep), 'k');
        plot(X_p(2:Np), Y_p(2:Np), 'b.', 'MarkerSize', 10)
        hold off
        
        set(gca,'Ylim',[min_y max_y])
        set(gca,'Xlim',[min_x max_x])
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer')
        frame = getframe(gcf); 
        if writevideo == true      
            writeVideo(writerObj,frame);
        end
    end
    
    timestep = timestep + 1;
end

%Close video writer
if writevideo == true
    close(writerObj);
end

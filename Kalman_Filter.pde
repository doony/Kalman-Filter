// Live Kalman filter demonstration
// Don Moloney z3361969 UNSW
// This is entirely my own work, with various bits of help from google, processing.org, youtube videos and online guides.
// This was pretty helpful: (Kalman Filter Basics, Monash Clayton) http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf
// And this: (Mean Vector and Covariance Matrix) http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm
// The "student dave" youtube and his ninja/quail code was good for refreshing the basic principles of the Kalman Filter. 

PFont f;  
float textBoxWidth = 0.4;            // How much of the screen to dedicate to text
float GNSS_sigma = 10.0;             // How noisy is the GNSS data? (stdev of measurement in meters)      (note that randomGaussian() returns mean 0 stdev 1)
float Accel_Noise = 10.0;            // How noisy is our accelerometer? 1.0 means a stdev of 1.0m/s/s for each measurement
float speed_limit = 50;              // Speed limit of our subject
float G_limit = 10.0;                // Maximum acceleration in m/s/s's of our subject - ability to start/stop/turn, each axis
float wobble = 1.0;                  // Multiplier for each change in direction (affects accel, won't exceed G_limit). Not related to Ez, Ex.
float sample_rateH = 20.0;           // State_act refresh rate, cycles per second                      (Time domain 20Hz)
float sample_rate  =  1.0;           // GNSS sample rate - cycles per second                           (Time domain  1Hz) 
float dtH = 1.0f / sample_rateH;     // Time interval between state_act update, seconds                (Time domain 20Hz)
float dt  = 1.0f / sample_rate;      // Time interval between GNSS updates, seconds.                   (Time domain  1Hz)
int custom_Q_R = 1;                  // Set this to one if you wish to define a custom Q and R, set below. Otherwise the values will be estimated based on data acquired over a period of covariance_sample seconds.
int covariance_sample = 7;           // How many seconds of data used to calibrate covariance
int debug_mode = 0;                  // Debug mode can be switched on/off with +/- keys
int error_graph_length = 100 ;       // Length of error graph
int background_fade = 2;             // How quickly the data fades from the main window. 0 means the data doesn't fade. 10 it fades almost immediately. 

float[] Q_custom = {Accel_Noise*Accel_Noise * dt*dt*dt*dt/4 , Accel_Noise * Accel_Noise *dt*dt, Accel_Noise*Accel_Noise * dt*dt*dt*dt/4 , Accel_Noise * Accel_Noise *dt*dt};
float[] R_custom = {GNSS_sigma * GNSS_sigma, GNSS_sigma * GNSS_sigma, GNSS_sigma * GNSS_sigma, GNSS_sigma * GNSS_sigma};
float[] state_act  = {450,0,450,0};         // Actual subject position {X, X_velocity, Y, Y_velocity}                                 (Time domain 20Hz)
float[] state_temp1 = new float[4];
float[] state_temp2 = new float[4];
float[] u_state_input_act =  {0,0};         // Acceleration (state input) vector x'' and y'' actual, used to generate boat path       (Time domain 20Hz)

float[] u_state_input_meas = {0,0};         // Acceleration (state input) vector x'' and y'' according to INS ? Noisy?                (Time domain  1Hz) 
float[] X_k_k      = {450,0,450,0};         // State Estimate(k|k)     {X, X_velocity, Y, Y_velocity}
float[] X_k_k_temp = {450,0,450,0};
float[] X_k1_k     = {450,0,450,0};         // State Estimate(k+1|k)   {X, X_velocity, Y, Y_velocity}
float[] X_k1_k1    = {450,0,450,0};         // State Estimate(k+1|k+1) {X, X_velocity, Y, Y_velocity}

float[][] F_state_trans_dtH = { {   1.0,   dtH,     0,     0 },        // State Transition Matrix F                                   (Time domain 20Hz) 
                                {     0,   1.0,     0,     0 },
                                {     0,     0,   1.0,   dtH },
                                {     0,     0,     0,   1.0 } };
                                
float[][] F_state_trans_dt =  { {   1.0,    dt,     0,     0 },        // State Transition Matrix F                                   (Time domain  1Hz) 
                                {     0,   1.0,     0,     0 },
                                {     0,     0,   1.0,    dt },
                                {     0,     0,     0,   1.0 } };
                          
float[][] G_input_control_dtH=  { {0.5*dtH*dtH,            0 },        // Input Control Matrix                                        (Time domain 20Hz) 
                                {          dtH,            0 },
                                {            0,  0.5*dtH*dtH },
                                {            0,          dtH } };
                               
float[][] G_input_control_dt =  { {  0.5*dt*dt,            0 },        // Input Control Matrix                                        (Time domain  1Hz) 
                                {           dt,            0 },
                                {            0,    0.5*dt*dt },
                                {            0,         dt } };       
                                
                                                 
float[][] H_meas_mx =         { {     1,      0,     0,     0 },        // Measurement Matrix - which elements to measure             (Time domain 1Hz) 
                                {     0,      1,     0,     0 },
                                {     0,      0,     1,     0 },
                                {     0,      0,     0,     1 }};
                            
float[] Z_k1                 =  {450,0,450,0};
float[] Z_k1_k               =  {450,0,450,0};

float[] V_k1 = {0,0,0,0};
float[] V_k1_temp = {0,0,0,0}; 
float[] Gu = new float[G_input_control_dt.length];
float[][] W_kalman_gain =            {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] P_k_k =                    {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] P_k1_k =                   {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] P_k1_k1 =                  {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] FPF_temp =                 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; // Used in Step 1 of State Covariance estimation
float[][] S_k1 =                     {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] HP_temp    =               {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; // Used in Step 2 of State Covariance estimation

float[][] X_k1_k_history = new float[covariance_sample][4];
float[]   X_k1_k_mean    =           {0,0,0,0};
float[][] Z_k1_history   = new float[covariance_sample][4];
float[]   Z_k1_mean      =           {0,0,0,0};
float[][] Q_covariance   =           {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
float[][] R_covariance   =           {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

float[] A_mag = {0};
float[] X_error = {0};
float[] X_v_error = {0};
float[] Y_error = {0};
float[] Y_v_error = {0};
float Pos_error_mean = 0;
float V_error_mean = 0;
float Pos_error_sd = 0;
float V_error_sd = 0;
 
int lastmillis;
float timesnap = second() + 60*minute() + 60*60*hour();
float lastsnap = second() + 60*minute() + 60*60*hour();
float timesnap0 = second() + 60*minute() + 60*60*hour();

public void setup() {
  size(1800,1000);               // Size of work area
  stroke(204,102,0);              // Line colour set to orange, redundant
                         // Anti-aliasing
  background(0);                  // Black background, redundant
f = createFont("Arial",16,true);  // Font selection for text area
if (custom_Q_R == 1) { covariance_sample =1;}
}

public void draw() {
  
  // Define keyboard controls
  if (keyPressed == true) {
    if (key == '1'){ GNSS_sigma = GNSS_sigma -0.1f; if (GNSS_sigma < 0){ GNSS_sigma = 0;} } 
    if (key == '2'){ GNSS_sigma = GNSS_sigma +0.1f;} 
    if (key == 'q'){ Accel_Noise /= 1.0233f; if (Accel_Noise < 0){ Accel_Noise = 0;} } 
    if (key == 'w'){ Accel_Noise *=  1.0233f;} 
    if (key == 'a'){ G_limit /=  1.0233f; if (G_limit < 0){ G_limit = 0;} } 
    if (key == 's'){ G_limit *= 1.0233f;} 
    if (key == 'z'){ speed_limit /= 1.0233f; if (speed_limit < 0){ speed_limit = 0;} } 
    if (key == 'x'){ speed_limit *=  1.0233f;} 
    if (key == '3'){ Q_covariance[0][0] /= 1.0233f;}
    if (key == '4'){ Q_covariance[0][0] *= 1.0233f;}
    if (key == 'e'){ Q_covariance[1][1] /= 1.0233f; }
    if (key == 'r'){ Q_covariance[1][1] *= 1.0233f; }
    if (key == 'd'){ Q_covariance[2][2] /= 1.0233f; }    
    if (key == 'f'){ Q_covariance[2][2] *= 1.0233f; }    
    if (key == 'c'){ Q_covariance[3][3] /= 1.0233f; }    
    if (key == 'v'){ Q_covariance[3][3] *= 1.0233f; }
  //  if (key == 'r'){ reset_all(); setup(); }
   if (key == '5'){ R_covariance[0][0] /= 1.0233f; }
   if (key == '6'){ R_covariance[0][0] *= 1.0233f; }
   if (key == 't'){ R_covariance[1][1] /= 1.0233f; }
   if (key == 'y'){ R_covariance[1][1] *= 1.0233f; }
   if (key == 'g'){ R_covariance[2][2] /= 1.0233f; }    
   if (key == 'h'){ R_covariance[2][2] *= 1.0233f; }    
   if (key == 'b'){ R_covariance[3][3] /= 1.0233f; }    
   if (key == 'n'){ R_covariance[3][3] *= 1.0233f; }
   if (key == '-'){ debug_mode=0;println(debug_mode);}
   if (key == '+'){ debug_mode=1;println(debug_mode);}
  }
    
  //Draw Layout, Background, grid lines
  stroke(204,102,0); fill(0,background_fade); rect(0,0,floor(width*(1-textBoxWidth)),height-1); //Main Window
  stroke(204,102,0); fill(0); rect(floor(width*(1-textBoxWidth)),0,width-1,height-1);  //Text Area
  stroke(80); for (int i = 1; i < height / 100; i++){line(0,100*i,width*(1-textBoxWidth),100*i);} for (int i = 1; i < width*(1-textBoxWidth) / 100; i++){  line(100*i,0,100*i,height); } // Grid Lines
  
  //Generate the boat's path - this occurs at higher frequency (sample_rateH = 20Hz) than our kalman filter is sampling (1 hz)
  
  noStroke();
  fill(255);
  
  for (int i = 0; i < state_temp1.length; i++) {
    state_temp1[i] = 0;
    state_temp2[i] = 0;
  }
  
  // State transition at 20Hz - iterating the boat state
  cross_mult(F_state_trans_dtH,state_act,state_temp1);
  cross_mult(G_input_control_dtH,u_state_input_act,state_temp2);
  
  
  for (int i=0;i<state_act.length; i++){
  state_act[i] = state_temp1[i] + state_temp2[i]; 
  }
  if (debug_mode == 1) {
  text("state_temp1",width-680,230);
  text("state_temp2",width-580,230);
  text("state_act",  width-480,230);
  for (int i=0; i < state_temp1.length; i++) {
     text(state_temp1[i],width-680,250+i*20);
     text(state_temp2[i],width-580,250+i*20);
     text(state_act[i],  width-480,250+i*20);
  } 
  }
  //Implementing a top speed for the vehicle
  if (abs(state_act[1]) > speed_limit) { state_act[1] *= speed_limit/abs(state_act[1]);}
  if (abs(state_act[3]) > speed_limit) { state_act[3] *= speed_limit/abs(state_act[3]);}
  
  //------------------START AUTOPILOT CODE (20Hz)------------------START AUTOPILOT CODE (20Hz)------------------START AUTOPILOT CODE (20Hz)------------------START AUTOPILOT CODE (20Hz)------------------START AUTOPILOT CODE (20Hz)
                                                                                                                 // The boat randomly decides to change direction. 
  float change_dir = random(10);
  if (timesnap == timesnap0+covariance_sample + 1) {
     u_state_input_act[0] = wobble * (random(2*G_limit)-G_limit); u_state_input_act[1] = wobble * (random(2*G_limit)-G_limit);
    if (abs(u_state_input_act[0]) > G_limit) { u_state_input_act[0] *= G_limit/abs(u_state_input_act[0]); }
    if (abs(u_state_input_act[1]) > G_limit) { u_state_input_act[1] *= G_limit/abs(u_state_input_act[1]); }
  }
  if (change_dir < 0.1f && timesnap > timesnap0 + covariance_sample) {
    u_state_input_act[0] = wobble * (random(2*G_limit)-G_limit); u_state_input_act[1] = wobble * (random(2*G_limit)-G_limit);
    if (abs(u_state_input_act[0]) > G_limit) { u_state_input_act[0] *= G_limit/abs(u_state_input_act[0]); }
    if (abs(u_state_input_act[1]) > G_limit) { u_state_input_act[1] *= G_limit/abs(u_state_input_act[1]); }
  }
                                                                                                                 // Bouncing off the walls
  if (state_act[0] < 1) { state_act[1] *= -1;} if (state_act[0] > width*(1-textBoxWidth)) { state_act[1] *=-1 ; state_act[0] = width*(1-textBoxWidth); }
  if (state_act[2] < 1) { state_act[3] *=-1;}  if (state_act[2] > height-1)               { state_act[3] *= -1; state_act[2] = height-1;               }
      
                                                                                                                 // When autopilot gets close to the walls it accelerates away (X)
  if (state_act[0] < 200 && state_act[1] < 0.0f) { u_state_input_act[0] = 0.5f*G_limit;                            // If within 200m of left wall, accelerate away at 1/2 throttle. 
  if (state_act[0] < 100 && state_act[1] < 0.0f) { u_state_input_act[0] = G_limit;    }}                          // If within 100m of left wall, accelerate away at 1.0 throttle. 
  if (state_act[0] > width*(1-textBoxWidth)-200 && state_act[1] > 0.0f) { u_state_input_act[0] = -0.5f*G_limit;    //ditto right wall
  if (state_act[0] > width*(1-textBoxWidth)-100 && state_act[1] > 0.0f) { u_state_input_act[0] = -G_limit;    }}  //ditto right wall
  if (state_act[2] < 200 && state_act[3] < 0.0f) { u_state_input_act[1] = 0.5f*G_limit;                            //ditto top wall
  if (state_act[2] < 100 && state_act[3] < 0.0f) { u_state_input_act[1] = G_limit;    }}                          //ditto top wall
  if (state_act[2] > height-200                 && state_act[3] > 0.0f) { u_state_input_act[1] = -0.5f*G_limit;    //ditto bottom wall
  if (state_act[2] > height-100                 && state_act[3] > 0.0f) { u_state_input_act[1] = -G_limit;    }}  //ditto bottom wall
  
  text("Auto pilot ON",5+floor(width*(1-textBoxWidth)),height-200);
  //------------------END AUTOPILOT CODE (20Hz)---------------------END AUTOPILOT CODE (20Hz)---------------------END AUTOPILOT CODE (20Hz)---------------------END AUTOPILOT CODE (20Hz)---------------------END AUTOPILOT CODE (20Hz)
  

  text("Velocity " + floor(sqrt(state_act[1]* state_act[1] + state_act[3]*state_act[3])),5+floor(width*(1-textBoxWidth)),height-160);
  if (state_act[1] == 0) {
    if (state_act[3] >0) {text("Heading 90",5+floor(width*(1-textBoxWidth)),height-180);}
    if (state_act[3] <=0) {text("Heading -90",5+floor(width*(1-textBoxWidth)),height-180);}
  }
    else {  
      if (state_act[1] < 0) {
        if (state_act[3] < 0) { text("Heading " + floor(-180+(1/3.14159f)*180*atan(state_act[3]/state_act[1])),5+floor(width*(1-textBoxWidth)),height-180);}
        if (state_act[3] >= 0) { text("Heading " + floor(180-(1/3.14159f)*180*atan(state_act[3]/state_act[1])),5+floor(width*(1-textBoxWidth)),height-180);}
        }
      if (state_act[1] > 0) {
        if (state_act[3] < 0) { text("Heading " + floor((1/3.14159f)*180*atan(state_act[3]/state_act[1])),5+floor(width*(1-textBoxWidth)),height-180);}
        if (state_act[3] >= 0) { text("Heading " + floor((1/3.14159f)*180*atan(state_act[3]/state_act[1])),5+floor(width*(1-textBoxWidth)),height-180);}
      }
    }
  
  
  
  //Velocity (green) and Accel (red) vector gauge
  fill(0); stroke(80); ellipse(floor(width*(1-textBoxWidth))+330,70,100,100);
  stroke(255,0,0); arrow(floor(width*(1-textBoxWidth))+330,70,floor(width*(1-textBoxWidth))+330+state_act[1]*30/speed_limit,70 + state_act[3]*30/speed_limit); // Velocity
  stroke(0,255,0); arrow(floor(width*(1-textBoxWidth))+330,70,floor(width*(1-textBoxWidth))+330+u_state_input_act[0]*30/G_limit,70 + u_state_input_act[1]*30/G_limit); // Accel
  fill(255,0,0); text("Velocity vector",floor(width*(1-textBoxWidth))+300,120);
  fill(0,255,0); text("Acceleration vector",floor(width*(1-textBoxWidth))+300,140);
    
  // Draw the boat's path! (white line)
  stroke(255);
  ellipse(state_act[0], state_act[2], 2, 2);

  // Start 1Hz time domain-----------------Start 1Hz time domain-----------------Start 1Hz time domain-----------------Start 1Hz time domain-----------------Start 1Hz time domain-----------------Start 1Hz time domain

  timesnap = second() + 60*minute() + 60*60*hour();
  if (timesnap == lastsnap + 1) {                               //START 1Hz DOMAIN  (sample_rate = 1Hz)
    lastsnap = second() + 60*minute() + 60*60*hour();
        
        
        
        //------------- BEGIN UPDATED STATE ESTIMATION/PREDICTION------------- BEGIN UPDATED STATE ESTIMATION/PREDICTION------------- BEGIN UPDATED STATE ESTIMATION/PREDICTION------------- BEGIN UPDATED STATE ESTIMATION/PREDICTION
    
    X_k_k = X_k1_k1;
    
                                                                                          // Start step 1 - Current Estimate X_k1_k
    u_state_input_meas[0] = u_state_input_act[0] + Accel_Noise * randomGaussian();        // Noisy accelerometer data
    u_state_input_meas[1] = u_state_input_act[1] + Accel_Noise * randomGaussian();        // Noisy accelerometer data
                                                                                          //
    if (abs(u_state_input_meas[0]) > G_limit) {                                           //
      u_state_input_meas[0] *= G_limit/abs(u_state_input_meas[0]);}                       // Enforcing acceleration limit on u, will affect X_k1_k
    if (abs(u_state_input_meas[1]) > G_limit) {                                           //
      u_state_input_meas[1] *= G_limit/abs(u_state_input_meas[1]);}                       //
                                                                                          //
    cross_mult(F_state_trans_dt,    X_k_k,    X_k_k_temp);                                // First term in the State Predict ( x = Fx + Gu) 
    cross_mult(G_input_control_dt,  u_state_input_meas, Gu  );                            // second term in the State Predict ( x = Fx + Gu)
                                                                                          //
    for (int i=0;i<X_k1_k.length;i++) {                                                   //
      X_k1_k[i] = X_k_k_temp[i] + Gu[i];                                                  // Finished Step 1 of State Estimation - we now have our initial Estimate (X_k1_k) 
    }
   if (timesnap < timesnap0 + covariance_sample + 1) {
      for (int i = 0; i < X_k1_k_history[0].length; i++){
        X_k1_k_history[floor(timesnap-timesnap0)-1][i] = X_k1_k[i];                       // We're filling up our history store
      }
    }
           
   if (abs(X_k1_k[1]) > speed_limit) { X_k1_k[1] *= speed_limit/abs(X_k1_k[1]);}          // Enforcing speed limit on X_k1_k
   if (abs(X_k1_k[3]) > speed_limit) { X_k1_k[3] *= speed_limit/abs(X_k1_k[3]);}
   
                                                                                          // Start Step 2 State Estimation... what is the predicted state? Looking for Z_k1_k
    cross_mult(H_meas_mx,X_k1_k, Z_k1_k);                                                 // Done! We have Z_k1_k our predicted state.
    noStroke(); fill(0,0,255); ellipse(Z_k1_k[0],Z_k1_k[2], 10,10);                       // Draw a Blue circle for Z_k1_k  (initial predicted state)
                                                                                          // Finished Step 2
    
                                                                               // Start Step 3 - what is the measurement residual? (V_k1)
                                                                               // First need the Z_k1 or GPS data..
    cross_mult(H_meas_mx,state_act,Z_k1);                                      // Recording GPS Data - this is the z(k+1) term or Z_k1 vector
    for (int i = 0; i < Z_k1.length; i ++) {                                   //
      Z_k1[i] += GNSS_sigma * randomGaussian() * H_meas_mx[i][i];              // GPS data reflects the current actual state + random gaussian distribution (adjustable via GNSS_sigma)
    }                                                                          //
    if (timesnap < timesnap0 + covariance_sample + 1) {                        //
      for (int i = 0; i < Z_k1_history[0].length; i++){                        //
        Z_k1_history[floor(timesnap-timesnap0)-1][i] = Z_k1[i];                // Fill our sample/historical data set up (Z_k1)
      }                                                                        //
    }                                                                          //
                                                                               // Draw the GPS data on graph! 
    noStroke(); fill(0,180,0); ellipse(Z_k1[0],Z_k1[2], 10,10);                // Green Circle: raw GNSS measurement
    stroke(255,0,0); fill(0,180,0);                                            //
    arrow(Z_k1[0],Z_k1[2],Z_k1[0]+Z_k1[1],Z_k1[2]+Z_k1[3]);                    // Red arrows: GNSS measured velocity vector
                                                                               //
                                                                               // Now calculate the residual V_k1
    for (int i = 0;i<V_k1.length; i++) {                                       //
      V_k1[i] = Z_k1[i] - Z_k1_k[i];                                           //
    }                                                                          // End Step 3 - we have measurement residual V_k1
    
                                                                               // Can't do Step 4 until we have done our state covariance estimates, which will provide us our Kalman Gain (W)
                                                                     
    
//-------------START COVARIANCE ESTIMATION-------------START COVARIANCE ESTIMATION-------------START COVARIANCE ESTIMATION-------------START COVARIANCE ESTIMATION-------------START COVARIANCE ESTIMATION-------------START COVARIANCE ESTIMATION
    
    for (int i = 0; i < P_k_k.length; i++) {                                               //
    for (int j = 0; j < P_k_k.length; j++) {                                               //
      P_k_k[i][j] = P_k1_k1[i][j];                                                         //
    }} 
    
    cross_mult_sq(F_state_trans_dt, P_k_k, FPF_temp);                                       // Part 1 COVARIANCE ESTIMATION:
    matrix_transpose(F_state_trans_dt);                                                     //  State Prediction Covariance (P_k1_k)
    cross_mult_sq(FPF_temp, F_state_trans_dt, FPF_temp);                                    //
    matrix_transpose(F_state_trans_dt);                                                     //
    for (int i = 0; i < P_k1_k.length; i++) {                                               //
    for (int j = 0; j < P_k1_k.length; j++) {                                               //
      P_k1_k[i][j] = FPF_temp[i][j] + Q_covariance[i][j];                                   //
    }}                                                                                      // Finished Part 1 - we have P_k1_k
    
    cross_mult_sq(H_meas_mx, P_k1_k, HP_temp);                                              // Part 2 COVARIANCE ESTIMATION:
    matrix_transpose(H_meas_mx);                                                            //  Measurement Prediction Covariance (S_k1)
    cross_mult_sq(HP_temp, H_meas_mx, HP_temp);                                             // 
    matrix_transpose(H_meas_mx);                                                            //
    for (int i = 0; i < S_k1.length; i++) {                                                 // 
    for (int j = 0; j < S_k1.length; j++) {                                                 // 
      S_k1[i][j] = HP_temp[i][j] + R_covariance[i][j];                                      // 
    }}                                                                                      // Finished Part 2 - we have S_k1
    
    matrix_transpose(H_meas_mx);                                                            // Part 3 COVARIANCE ESTIMATION: 
    cross_mult_sq(P_k1_k, H_meas_mx, HP_temp);                                              //  Calculating Filter Gain (W)
    matrix_transpose(H_meas_mx);                                                            //
    matrix_inv(S_k1);                                                                       //
    cross_mult_sq(HP_temp, S_k1, W_kalman_gain);                                            //
    matrix_inv(S_k1);                                                                       // Finished Part 3 - we have W_kalman_gain
    
    cross_mult_sq(W_kalman_gain, S_k1, HP_temp);                                            // Part 4 COVARIANCE ESTIMATION: 
    matrix_transpose(W_kalman_gain);                                                        //  Calculating Updated State Covariance (P_k1_k1)
    cross_mult_sq(HP_temp,W_kalman_gain,HP_temp);                                           // 
    matrix_transpose(W_kalman_gain);                                                        // 
    for (int i = 0; i < P_k1_k1.length; i++) {                                              // 
    for (int j = 0; j < P_k1_k1.length; j++) {                                              // 
      P_k1_k1[i][j] = P_k1_k[i][j] - HP_temp[i][j];                                         // Part 4 complete, we have P_k1_k1 - updated state covariance. 
    }}

//-------------END COVARIANCE ESTIMATION-------------END COVARIANCE ESTIMATION-------------END COVARIANCE ESTIMATION-------------END COVARIANCE ESTIMATION-------------END COVARIANCE ESTIMATION-------------END COVARIANCE ESTIMATION
    
    
    cross_mult(W_kalman_gain,V_k1,V_k1_temp);                                                           // Step 4 - Calculate the updated estimate. X_k1_k = X_k1_k + W_kalman_gain * (Z_k1 - H_meas_mx * X_k1_k) 
    for (int i = 0; i<X_k1_k1.length;i++) {                                                             //
      X_k1_k1[i] = X_k1_k[i] + V_k1_temp[i];                                                            // Step 4 complete. We have our updated state estimate X_k1_k1
    }
    stroke(0,255,255); fill(0,255,255);     ellipse(X_k1_k1[0],X_k1_k1[2], 5,5);                        // Plotting the updated state estimate. We'll use a little turquoise circle for X_k1_k1
    
    if (X_error.length < error_graph_length){ X_error = append(X_error,X_k1_k1[0]-state_act[0]);}       // This is calculating the error in the updated state estimate and recording it to a ~100 long rolling vector
      else { shift_vector(X_error); X_error[error_graph_length-1] = X_k1_k1[0]-state_act[0];}
    if (X_v_error.length < error_graph_length){ X_v_error= append(X_v_error,X_k1_k1[1]-state_act[1]);} 
      else { shift_vector(X_v_error); X_v_error[error_graph_length-1] = X_k1_k1[1]-state_act[1];}
    if (Y_error.length < error_graph_length){ Y_error = append(Y_error,X_k1_k1[2]-state_act[2]);} 
      else { shift_vector(Y_error); Y_error[error_graph_length-1] = X_k1_k1[2]-state_act[2];}
    if (Y_v_error.length < error_graph_length){ Y_v_error = append(Y_v_error,X_k1_k1[3]-state_act[3]); } 
      else { shift_vector(Y_v_error); Y_v_error[error_graph_length-1] = X_k1_k1[3]-state_act[3];}
    if (A_mag.length < error_graph_length) { A_mag = append(A_mag,sqrt(u_state_input_act[0]*u_state_input_act[0]+u_state_input_act[1]*u_state_input_act[1]));} 
      else { shift_vector(A_mag); A_mag[error_graph_length-1] = sqrt(u_state_input_act[0]*u_state_input_act[0]+u_state_input_act[1]*u_state_input_act[1]);}
    
   
     if (timesnap == timesnap0 + covariance_sample) {                                                    // This code is for automatically generating Q and R diagonal matrices based on initial variance. 
       for (int i = 0; i<X_k1_k_history[0].length; i++) {                                                // That is, the subject remains stationary for covariance_sample seconds
         for (int j = 0; j<X_k1_k_history.length; j++) {                                                 // .. and calculates the elements of Q,R based on the variance of the INS and GNSS data received.
           X_k1_k_mean[i] += X_k1_k_history[j][i];                                                       //
           Z_k1_mean[i]     += Z_k1_history[j][i];                                                       //
         }                                                                                               //
         X_k1_k_mean[i] /= covariance_sample;                                                            // Finding our state mean values
         Z_k1_mean[i] /= covariance_sample;                                                              //
         }                                                                                               //
         for (int i=0;i<Q_covariance.length;i++) {                                                       // Calculate the variances for each state element e.g. X, X_velocity etc (i,i - along the diagonal)
           for (int j = 0; j <X_k1_k_history.length;j++){                                                //
             Q_covariance[i][i] += (X_k1_k_history[j][i]-X_k1_k_mean[i])*(X_k1_k_history[j][i]-X_k1_k_mean[i]); 
             R_covariance[i][i] += (  Z_k1_history[j][i]-  Z_k1_mean[i])*(  Z_k1_history[j][i]-  Z_k1_mean[i]);
           }                                                                                             //
         Q_covariance[i][i] /= (X_k1_k_history.length-1);                                                //
         R_covariance[i][i] /= (  Z_k1_history.length-1);                                                //
      if (custom_Q_R == 1) {                                                                             // However we have the option instead to set custom_Q_R to 1 and defining our Q,R matrices in the setup above.
        Q_covariance[i][i] = Q_custom[i];                                                                // .. or adjusting Q,R on the fly within the program using keyboard controls 
        R_covariance[i][i] = R_custom[i];                                                                // End of Covariance sampling code
      }
  
    }
      
        //calculate the covariance (CORRELATIONS) between variables at i/j - i.e. not on the diagonals  - This seems to spaz out the predictions so it's been commented out. 
        //for (int i=0;i<Q_covariance.length-1;i++) {
        // for (int j = i+1; j <Q_covariance.length;j++) {
        //   for (int k = 0; k<X_k1_k_history.length; k++) {
        //   Q_covariance[i][j] += (X_k1_k_history[k][i]-X_k1_k_mean[i])*(X_k1_k_history[k][j]-X_k1_k_mean[j])/(X_k1_k_history.length-1);
        //   R_covariance[i][j] += (  Z_k1_history[k][i]-  Z_k1_mean[i])*(  Z_k1_history[k][j]-  Z_k1_mean[j])/(  Z_k1_history.length-1);
        //   Q_covariance[j][i] += (X_k1_k_history[k][i]-X_k1_k_mean[i])*(X_k1_k_history[k][j]-X_k1_k_mean[j])/(X_k1_k_history.length-1);
        //   R_covariance[j][i] += (  Z_k1_history[k][i]-  Z_k1_mean[i])*(  Z_k1_history[k][j]-  Z_k1_mean[j])/(  Z_k1_history.length-1);
        //   }
        //   Q_covariance[j][i] = Q_covariance[i][j];
        //   R_covariance[j][i] = R_covariance[i][j];
            
        // }
        //}
        timesnap = 0;
        }                                                                                        // End calculation of Q and R
    

    
  }        // -------------------End 1Hz time domain-------------------End 1Hz time domain-------------------End 1Hz time domain-------------------End 1Hz time domain
  
  
  if (debug_mode == 0) {               //Analysis view items - displays graphs of error and shows Q,R controls
    stroke(255);
    line(width*(1-textBoxWidth),height-500,width*(1-textBoxWidth)+400,height-500);stroke(100);
    line(width*(1-textBoxWidth),height-540,width*(1-textBoxWidth)+400,height-540);
    line(width*(1-textBoxWidth),height-580,width*(1-textBoxWidth)+400,height-580);
    line(width*(1-textBoxWidth),height-620,width*(1-textBoxWidth)+400,height-620);
    line(width*(1-textBoxWidth),height-660,width*(1-textBoxWidth)+400,height-660);stroke(255);
    //for (int i = 1; i<X_error.length; i++){             // Plot seperate axis - proved to be too cluttered (4 data sets on one graph) - switched to position/velocity only
    //  fill(255,50,0); stroke(255,50,0);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(X_error[i]),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(X_error[i-1]));
    ////ellipse(width-200,height-500,10,10+X_error[0]);
    //}
    //for (int i = 1; i<X_v_error.length; i++){
    //  fill(0,255,50); stroke(0,255,50);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(X_v_error[i]),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(X_v_error[i-1]));
    ////ellipse(width-200,height-500,10,10+X_error[0]);
    //}
    //for (int i = 1; i<Y_error.length; i++){
    //  fill(255,100,0); stroke(255,100,0);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(Y_error[i]),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(Y_error[i-1]));
    ////ellipse(width-200,height-500,10,10+X_error[0]);
    //}
    //for (int i = 1; i<Y_v_error.length; i++){
    //  fill(0,255,100); stroke(0,255,200);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(Y_v_error[i]),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(Y_v_error[i-1]));
    ////ellipse(width-200,height-500,10,10+X_error[0]);
    //}
    Pos_error_mean = 0;        // Mean position error(magnitude)
    Pos_error_sd = 0;          // Standard deviation of position error (error magnitude.. not SD from the 'actual')
    V_error_mean = 0;          // Velocity error mean
    V_error_sd = 0;            // velocity error sd
    
     for (int i = 1; i<X_v_error.length; i++){
      fill(255,50,0); stroke(255,50,0);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(sqrt(X_v_error[i]*X_v_error[i]+Y_v_error[i]*Y_v_error[i])),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(sqrt(X_v_error[i-1]*X_v_error[i-1]+Y_v_error[i-1]*Y_v_error[i-1])));
      V_error_mean += sqrt(X_v_error[i]*X_v_error[i]+Y_v_error[i]*Y_v_error[i]);
    //ellipse(width-200,height-500,10,10+X_error[0]);
    }
    
    for (int i = 1; i<X_error.length; i++){
      fill(0,255,50); stroke(0,255,50);  line(width*(1-textBoxWidth)+4*i,height-500-5*abs(sqrt(X_error[i]*X_error[i]+Y_error[i]*Y_error[i])),width*(1-textBoxWidth)+4*(i-1),height-500-5*abs(sqrt(X_error[i-1]*X_error[i-1]+Y_error[i-1]*Y_error[i-1])));
      Pos_error_mean += sqrt(X_error[i]*X_error[i]+Y_error[i]*Y_error[i]);
    //ellipse(width-200,height-500,10,10+X_error[0]);
    }
    
    Pos_error_mean /= X_error.length;
    V_error_mean /= X_v_error.length;
    
    for (int i = 1; i<X_error.length; i++){
      Pos_error_sd += (sqrt(X_error[i]*X_error[i]+Y_error[i]*Y_error[i])-Pos_error_mean)*(sqrt(X_error[i]*X_error[i]+Y_error[i]*Y_error[i])-Pos_error_mean);
      V_error_sd +=   (sqrt(X_v_error[i]*X_v_error[i]+Y_v_error[i]*Y_v_error[i])-V_error_mean)*(sqrt(X_v_error[i]*X_v_error[i]+Y_v_error[i]*Y_v_error[i])-V_error_mean);
    }
    Pos_error_sd /= X_error.length; Pos_error_sd = sqrt(Pos_error_sd);
    V_error_sd /= X_v_error.length; V_error_sd = sqrt(V_error_sd);
    
    fill(255,0,0);
    text("Velocity Error",width*(1-textBoxWidth),height-660);
    text("mean="+V_error_mean+", stdv="+V_error_sd,width*(1-textBoxWidth)+150,height-660); fill(0,255,0);
    stroke(255,0,0); line(width*(1-textBoxWidth),height-500-5*V_error_mean,width*(1-textBoxWidth)+400,height-500-5*V_error_mean);
    text("Position Error",width*(1-textBoxWidth),height-680);
    text("mean="+Pos_error_mean+", stdv="+Pos_error_sd,width*(1-textBoxWidth)+150,height-680);
    stroke(0,255,0); line(width*(1-textBoxWidth),height-500-5*Pos_error_mean,width*(1-textBoxWidth)+400,height-500-5*Pos_error_mean);
    
    
    // ERROR vs ACCELERATION plot - to give an idea of relationship if any.
    text("Acceleration Magnitude",width-180,height-480);
    text("Position Error",        width-280,  height-620);fill(255,0,0);
    text("Velocity Error",        width-280,  height-600);fill(255);
    arrow(width-300,height-500,width-5,height-500);
    arrow(width-300,height-500, width-300,height-700);
    for (int i = 1; i<A_mag.length; i++){
      fill(0,255,0); stroke(0,255,0);  ellipse(width-300+300/1.5f/G_limit*abs(A_mag[i]),height-500-5*abs(sqrt(X_error[i]*X_error[i]+Y_error[i]*Y_error[i])),5,5);
      fill(255,0,0); stroke(255,0,0);  ellipse(width-300+300/1.5f/G_limit*abs(A_mag[i]),height-500-5*abs(sqrt(X_v_error[i]*X_v_error[i]+Y_v_error[i]*Y_v_error[i])),5,5);
    //ellipse(width-200,height-500,10,10+X_error[0]);
    }
    
  }
  
  
  
  
  
  
  if (debug_mode==1){                  // Displaying various matrices in debug mode helps me to understand when or why things are going wrong.
    
  fill(255);text("X_k1_k_history",width-200,230);fill(180);  
  for (int i=0; i < X_k1_k_history.length; i++) {
      for (int j=0; j < X_k1_k_history[0].length; j++) {
        text(floor(X_k1_k_history[i][j]),width-200+j*40,250+i*20);
  }}
  fill(255);
  text("Z_k1_history",width-200,570);fill(180);
  for (int i=0; i < Z_k1_history.length; i++) {
      for (int j=0; j < Z_k1_history[0].length; j++) {
        text(floor(Z_k1_history[i][j]),width-200+j*40,590+i*20);
  }}
  
  fill(255); text("X_k_k",width-680,330); fill(180);
  for (int i=0; i < X_k_k.length; i++) {
  text(X_k_k[i],width-680,350+i*20);
  }
  fill(255); text("X_k1_k",width-580,330); fill(180);
  for (int i=0; i < X_k1_k.length; i++) {
    text(X_k1_k[i],width-580,350+i*20);
  }
  fill(255); text("Z_k1_k",width-480,330); fill(180);
  for (int i=0; i < Z_k1_k.length; i++) {
    text(Z_k1_k[i],width-480,350+i*20);
  }
  fill(255); text("V_k1",width-680,430); fill(180);
  for (int i=0; i < Z_k1_k.length; i++) {
    text(V_k1[i],width-680,450+i*20);
  }
  fill(255); text("P_k1_k1",width-580,430); fill(180);
  for (int i=0; i < Z_k1_k.length; i++) {
    text(P_k1_k1[i][i],width-580,450+i*20);
  }
  fill(255); text("W_kalman_gain",width-680,530); fill(180);
  for (int i=0; i < Z_k1_k.length; i++) {
    text(W_kalman_gain[i][i],width-680,550+i*20);
  }
  fill(255); text("X_k1_k1",width-480,430); fill(180);
  for (int i=0; i < X_k1_k1.length; i++) {
    text(X_k1_k1[i],width-480,450+i*20);
  }
  fill(255); text("H_meas_mx",width-380,330); fill(180);
  for (int i=0; i < H_meas_mx.length; i++) {
    text(H_meas_mx[i][i],width-380,350+i*20);
  }
  fill(255); text("u_state_input_meas",width-380,230); fill(180);
  for (int i=0; i < u_state_input_meas.length; i++) {
    text("Meas:" +u_state_input_meas[i],width-380,250+i*40);
    text("    Act:" + u_state_input_act[i],width-380,250+i*40+20);
  }   
  
  } //end debug mode
  
  
  //Q_Covariance
  fill(255);       text("Q_Covariance:",width*(1-textBoxWidth)+20,540+150*debug_mode);
  fill(280,255,0); text("(3) <                > (4)",width*(1-textBoxWidth)+20,560+150*debug_mode);
  fill(280,255,0); text("(E) <                > (R)",width*(1-textBoxWidth)+20,580+150*debug_mode);
  fill(280,255,0); text("(D) <                > (F)",width*(1-textBoxWidth)+20,600+150*debug_mode);
  fill(280,255,0); text("(C) <                > (V)",width*(1-textBoxWidth)+20,620+150*debug_mode);
  fill(255);
  for (int i=0; i < Q_covariance.length; i++) {
    text(Q_covariance[i][i],width*(1-textBoxWidth)+60,560+i*20+150*debug_mode);
  }
  
  //R_Covariance
  fill(255);       text("R_Covariance:",width*(1-textBoxWidth)+220,540+150*debug_mode);
  fill(280,255,0); text("(5) <                > (6)",width*(1-textBoxWidth)+220,560+150*debug_mode);
  fill(280,255,0); text("(T) <                > (Y)",width*(1-textBoxWidth)+220,580+150*debug_mode);
  fill(280,255,0); text("(G) <                > (H)",width*(1-textBoxWidth)+220,600+150*debug_mode);
  fill(280,255,0); text("(B) <                > (N)",width*(1-textBoxWidth)+220,620+150*debug_mode);
  fill(255);
  for (int i=0; i < R_covariance.length; i++) {
    text(R_covariance[i][i],width*(1-textBoxWidth)+260,560+i*20+150*debug_mode);
    
    
  }
  //draw data text at top of text area
  fill(280,255,0); text("Keyboard controls in YELLOW",width-300,40);
  fill(280,255,0); text("Debug Mode [+]",width-300,60);
  fill(280,255,0); text("Analysis Mode [-]",width-300,80);
  fill(255);
  textFont(f,16);
  text("GNSS Noise =          " + 0.1f*floor(GNSS_sigma*10)+" m",5+floor(width*(1-textBoxWidth)),20);
  text("INS  Noise =              " + 0.1f*floor(Accel_Noise*10)+" m/s2",5+floor(width*(1-textBoxWidth)),40);
  fill(280,255,0); text("  (1)<                  >(2)",5+floor(width*(1-textBoxWidth))+110,20);
  fill(280,255,0); text("  (Q)<                 >(W)",5+floor(width*(1-textBoxWidth))+110,40);
  //G_limit
  fill(255);  text("Accel Limit =              " + G_limit + " m/s2",5+floor(width*(1-textBoxWidth)),60);
  //Speed_limit
  fill(255);  text("Speed Limit  =           " + speed_limit + " m/s",5+floor(width*(1-textBoxWidth)),80);
  fill(280,255,0); text("  (A)<                  >(S)",5+floor(width*(1-textBoxWidth))+110,60);
  fill(280,255,0); text("  (Z)<                  >(X)",5+floor(width*(1-textBoxWidth))+110,80);
  //Actual/Estimate/Delta
  fill(180);
  text("        Actual  Est       Delta",5+floor(width*(1-textBoxWidth)),100); fill(255); 
  text("X        " +  floor(state_act[0]), 5+floor(width*(1-textBoxWidth)),  120); text( floor(X_k1_k1[0]),5+floor(width*(1-textBoxWidth))+90,  120); text(floor(X_k1_k1[0]-state_act[0]),5+floor(width*(1-textBoxWidth))+140,  120);
  text("Y        " +  floor(state_act[2]), 5+floor(width*(1-textBoxWidth)),  140); text( floor(X_k1_k1[2]),5+floor(width*(1-textBoxWidth))+90,  140); text(floor(X_k1_k1[2]-state_act[2]),5+floor(width*(1-textBoxWidth))+140,  140);
  text("Vel_X    " +  floor(state_act[1]), 5+floor(width*(1-textBoxWidth)),  160); text( floor(X_k1_k1[1]),5+floor(width*(1-textBoxWidth))+90,  160); text(floor(X_k1_k1[1]-state_act[1]),5+floor(width*(1-textBoxWidth))+140,  160);
  text("Vel_Y    " +  floor(state_act[3]), 5+floor(width*(1-textBoxWidth)),  180); text( floor(X_k1_k1[3]),5+floor(width*(1-textBoxWidth))+90,  180); text(floor(X_k1_k1[3]-state_act[3]),5+floor(width*(1-textBoxWidth))+140,  180);
  text("Acc_X    " +  u_state_input_act[0], 5+floor(width*(1-textBoxWidth)),  200);
  text("Acc_Y    " +  u_state_input_act[1], 5+floor(width*(1-textBoxWidth)),  220);

  
  //Draw LEGEND at bottom of text area
  text("State Estimate:",5+floor(width*(1-textBoxWidth)),height-80);   stroke(0,255,255); fill(0,255,255); ellipse(50+floor(width*(1-textBoxWidth))+100,height-85,5,5);
  fill(255);
  text("Actual:",  5+floor(width*(1-textBoxWidth)),height-60);        noStroke(); fill(255);     for (int i = 0; i < 20; i=i + 4) { ellipse(50+i+floor(width*(1-textBoxWidth))+ 96,height-65, 2, 2); }
  text("GNSS reading:",5+floor(width*(1-textBoxWidth)),height-40);   noStroke(); fill(0,180,0);   ellipse(50+floor(width*(1-textBoxWidth))+100,height-45,10,10); stroke(255,0,0);  arrow(50+floor(width*(1-textBoxWidth))+100,height-45,floor(width*(1-textBoxWidth))+170,height-45);
  fill(255);
  text("INS Predict:",5+floor(width*(1-textBoxWidth)),height-20);   noStroke(); fill(0,0,255);   ellipse(50+floor(width*(1-textBoxWidth))+100,height-25,10,10);
  
  // Program refresh set to 1rpt/50ms = 20Hz. Maximum refresh around 60Hz due to program time to repeat. println prints how many milliseconds between cycles. should be 50ms.
  while (millis() < lastmillis + 1000/sample_rateH){} println("[" + (millis()-lastmillis) + "]ms cycle should match [" + floor(dtH*1000) + "]ms period, since sample rate is set to ["+sample_rateH + "]Hz"); lastmillis = millis() ;
}



// Various functions used

public void cross_mult_n(float[][] A, float[] Q, float[] B) {    //Function multiplies A[NxN] by Q[Nx1]. N or 'int n' is the state vector length. Writes to input matrix B. (not really used - cross_mult works better)
  for (int i = 0; i < A.length; i++) {for (int j = 0; j < A.length; j++) {if (i==0 && j==0){B[0] =0; B[1]=0;} B[i] = B[i] + A[i][j]*Q[j];}}
}

public void cross_mult(float[][] A, float[] Q, float[] B){       //Function multiplies A[NxM] by Q[Mx1]. 
  int n = A.length; int m = A[0].length;
  for (int i = 0; i < n; i++) { B[i] = 0;} for (int i = 0; i < n; i++) { for (int j = 0; j < m; j++) { B[i] = B[i] + A[i][j]*Q[j];}}
}

public void cross_mult_sq(float[][] A, float[][] Q, float[][] B){       //Function multiplies two square matrices. 
int n = A.length;
  for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) { B[i][j] = 0; } }
  for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) { for (int k = 0; k < n; k++) { B[i][j] += A[i][k]*Q[k][j];}}}
}

public void matrix_inv(float[][] A){                             //Replaces 4x4 matrix with its inverse. See inv_matrix4x4_Working for comments. 
  float[] boo = new float[A.length * A.length]; 
  for (int i=0; i<A.length;i++) { for (int j=0; j<A.length;j++) { boo[A.length*i + j] = A[i][j]; } }
  PMatrix3D testmatrix = new PMatrix3D(); testmatrix.set(boo); testmatrix.invert(); testmatrix.get(boo); 
  for (int i=0; i < A.length;i++) { for (int j=0; j < A.length;j++) { A[i][j] = boo[A.length*i+j];  } }
}

public void matrix_transpose(float[][] A){                       //Replaces 4x4 matrix with its transposed matrix. 
  float[] boo = new float[A.length * A.length]; 
  for (int i=0; i<A.length;i++) { for (int j=0; j<A.length;j++) { boo[A.length*i + j] = A[i][j]; } }
  PMatrix3D testmatrix = new PMatrix3D(); testmatrix.set(boo); testmatrix.transpose(); testmatrix.get(boo);   
  for (int i=0; i < A.length;i++) { for (int j=0; j < A.length;j++) { A[i][j] = boo[A.length*i+j];  } }
}

public void arrow(float x1, float y1, float x2, float y2) {      //Draws an arrow from [1] to [2] - adapted from from user joridegoede https://processing.org/discourse/beta/num_1219607845.html
  line(x1, y1, x2, y2); pushMatrix(); translate(x2, y2); float a = atan2(x1-x2, y2-y1); rotate(a); line(0, 0, -5, -10); line(0, 0, 5, -10); popMatrix();
} 

public void reset_all() {
  for(int i = 0; i < state_act.length; i++){
    if (i==1 || i==3) {
      state_act[i]=0;X_k_k[i]=0;X_k_k_temp[i]=0;X_k1_k[i]=0;X_k1_k1[i]=0;Z_k1[i]=0;Z_k1_k[i]=0;V_k1[i]=0;V_k1_temp[i]=0;
    }
    if (i==0||i==2){
      state_act[i]=450;X_k_k[i]=450;X_k_k_temp[i]=450;X_k1_k[i]=450;X_k1_k1[i]=450;Z_k1[i]=450;Z_k1_k[i]=450;V_k1[i]=450;V_k1_temp[i]=450;
    }
    
    float timesnap = second() + 60*minute() + 60*60*hour();
    float lastsnap = second() + 60*minute() + 60*60*hour();
    float timesnap0 = second() + 60*minute() + 60*60*hour();
}
  stroke(204,102,0); fill(0,5); rect(0,0,floor(width*(1-textBoxWidth)),height-1); //Main Window
  stroke(204,102,0); fill(0); rect(floor(width*(1-textBoxWidth)),0,width-1,height-1);  //Text Area
  stroke(80); for (int i = 1; i < height / 100; i++){line(0,100*i,width*(1-textBoxWidth),100*i);} for (int i = 1; i < width*(1-textBoxWidth) / 100; i++){  line(100*i,0,100*i,height); } // Grid Lines
}

public void shift_vector(float[] A) {
 for (int i = 0;i<A.length-1;i++){
  A[i] = A[i+1]; 
 }}

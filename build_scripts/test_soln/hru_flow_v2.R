## TVDish solution
rm(list=ls())
graphics.off()

## inputs and time step
n <- 1100

q_sf_in <- rep(0,n)
q_sf_in[200:500] <- 0
q_sz_in <- rep(0.00001,n)
#q_sz_in[800:1000] <- 1
p <- rep(0.2/1000,n)
pet <- rep(0,n)


Dt <- 900
Dx <- 100
area <- 1000
t_d <- 360
s_rz_max <- 0.1*area
s_rz_0 <- 0.75

nstep <- 10
epsilon <- 0.001#1e-10 #0.001 #1e-10
q_uz_sz_init <- 0/area

## function to give outflow from the surface zone for an area
fqsf <- function(a,qin){
    max(0, 2*(1.1 * (a^(5/3))) - qin)
}
fasf <- function(qout,qin){ q <- (qout+qin)/2; (q/1.1)^(3/5) }

## function to give outflow from the surface zone for an area
fqsz <- function(a){ 1.1*exp(-0.2*a) }
fasz <- function(q){ -log(q/1.1)/0.2  }

s_sf <- s_rz <- s_uz <- s_sz <- q_sf <- q_sz <- rep(NA,n)
v_sf_rz<- v_rz_uz <- v_uz_sz <- rep(NA,n)
v_sf_in <- v_sf_out <- v_sz_in <- v_sz_out <- vp <- vep <- rep(NA,n)
q_sf_out <- q_sz_out <- rep(NA,n)

q_sz_max <- fqsz(0)


## initialise at steady state
qq = min( q_sz_in[1] , q_sz_max)
q_sf_in[1] <- q_sf_in[1] + q_sz_in[1] - qq
q_sz_in[1] <- qq

r_rz_uz <- r_sf_rz <- q_sf_in[1]
r_inj <- area * q_uz_sz_init

r_uz_sz <- min(r_rz_uz + r_inj, area/t_d)
r_uz_sz <- min(r_uz_sz, q_sz_max - q_sz_in[1])
q_sz_out[1] <- q_sz_in[1] + r_uz_sz

qq <- (q_sz_out[1] + q_sz_in[1])/2 ## since kinematic
s_sz[1] <- fasz(qq)*Dx
s_uz[1] = t_d * r_uz_sz * s_sz[1] / area
r_rz_uz = r_uz_sz - r_inj;

## solve rz
if( (r_sf_rz > 0.0) | (r_rz_uz < 0.0)  ){
    s_rz[1] = s_rz_max
}else{
    s_rz[1] = s_rz_max * s_rz_0
}

## balance flux through root zone
r_sf_rz = min( r_sf_rz , r_rz_uz )
q_sf_out[1] <- q_sf_in[1] - r_sf_rz
s_sf[1] <- Dx * fasf(q_sf_out[1],q_sf_in[1])

head(s_sf)
head(s_rz)
head(s_uz)
head(s_sz)
head(q_sz)
head(q_sf)
##head(v_sf_rz)
##head(v_rz_uz)
##head(v_uz_sz)

## loop
for(tt in 2:n){
    print(tt)
    qq = min( q_sz_in[tt] , q_sz_max)
    q_sf_in[tt] <- q_sf_in[tt] + q_sz_in[tt] - qq
    q_sz_in[tt] <- qq

    v_sf_in[tt] <- Dt*q_sf_in[tt]
    v_sz_in[tt] <- Dt*q_sz_in[tt]
    vp[tt] <- Dt*p[tt]*area
    vep[tt] <- 0

    vsfStep <- v_sf_in[tt]/nstep
    vszStep <- v_sz_in[tt]/nstep
    pStep <- vp[tt]/nstep
    DtStep <- Dt/nstep
    epStep <- pet[tt]*DtStep*area

    sz <- s_sz[tt-1]
    uz <- s_uz[tt-1]
    rz <- s_rz[tt-1]
    sf <- s_sf[tt-1]
    
    v_sf_rz[tt] <- v_rz_uz[tt] <- v_uz_sz[tt] <- 0
    v_sf_out[tt] <- v_sz_out[tt] <- 0
    
    for(ii in 1:nstep){
        ## downward pass
        vsfrz <- sf + vsfStep
        vrzuz <- max(0, rz + pStep - epStep + vsfrz - s_rz_max)
        ## solve sz
        ## test for saturation
        z <- 0
        vuzsz <- area*DtStep*min( (uz+vrzuz)/(t_d*z + area*DtStep), 1/t_d )
        qq <- min(q_sz_max,max(0,2*fqsz(z/Dx)-q_sz_in[tt]))
        Hz <- z - sz + vszStep + vuzsz - DtStep*qq
        #browser()
        if(Hz < 0){
            ##if(tt==200){ browser() }
            ## then solve
            rng <- c(z,sz+3*epsilon)
            ## expand
            flg <- TRUE
            ##print("expand")
            while(flg){
                z <- rng[2]
                vuzsz <- area*DtStep*min( (uz+vrzuz)/(t_d*z + area*DtStep), 1/t_d )
                qq <- min(q_sz_max,max(0,2*fqsz(z/Dx)-q_sz_in[tt]))
                Hz <- z - sz + vszStep + vuzsz - DtStep*qq
                if( Hz<0 ){
                    rng[1] <- rng[2]
                    rng[2] <- 2*rng[2]
                }else{
                    flg <- FALSE
                }
            }
            ## shrink
            ##print("shrink")
            while( diff(rng) > epsilon ){
                z <- (rng[1]+rng[2])/2
                vuzsz <- area*DtStep*min( (uz+vrzuz)/(t_d*z + area*DtStep), 1/t_d )
                qq <- min(q_sz_max,max(0,2*fqsz(z/Dx)-q_sz_in[tt]))
                Hz <- z - sz + vszStep + vuzsz - DtStep*qq
                if( Hz<= 0 ){rng[1] <- z}else{ rng[2] <- z}
            }
            z <- rng[2]
            print(rng)
        }

        ## upward pass
        #browser()
        qsz <- min(q_sz_max,max(0,2*fqsz(z/Dx)-q_sz_in[tt]))
        vuzsz <- sz - vszStep + DtStep*qsz - z
        #browser()
        sz <- z
        #browser()
        z <- min(sz, uz+vrzuz-vuzsz)
        vrzuz <- z + vuzsz - uz
        uz <- z
        ##vrzuz <- min( vrzuz, uz - vuzsz - sz )
        ##uz <- uz + vrzuz - vuzsz
        
        vsfrz <- min( vsfrz, s_rz_max - rz - pStep + epStep + vrzuz)
        rz <- (s_rz_max / (s_rz_max + epStep)) * (rz + pStep + vsfrz - vrzuz)

        ## surface
        rng <- c(0, sf + vsfStep - vsfrz)
        print("sf")
        while(diff(rng)>epsilon){
            w <- (rng[1]+rng[2])/2
            qq <- fqsf(w/Dx,q_sf_in[tt])
            Sw <- sf + vsfStep - vsfrz - (DtStep*qq) - w
            if(Sw<=0){rng[2] <- w}else{rng[1] <- w}
        }
        w <- rng[1]
        qsf <- (sf + vsfStep - vsfrz - w)/DtStep
        sf <- w

        v_sf_out[tt] <- v_sf_out[tt] + DtStep*qsf
        v_sz_out[tt] <- v_sz_out[tt] + DtStep*qsz
        v_sf_rz[tt] <- v_sf_rz[tt] <- vsfrz
        v_rz_uz[tt] <- v_rz_uz[tt] <- vrzuz
        v_uz_sz[tt] <- v_uz_sz[tt] <- vuzsz
    }

    s_sf[tt] <- sf
    s_rz[tt] <- rz
    s_uz[tt] <- uz
    s_sz[tt] <- sz
    q_sf_out[tt] <- qsf
    q_sz_out[tt] <- qsz
    
}

s <- s_sf + s_rz + s_uz - s_sz
vin <- v_sf_in + v_sz_in + vp
vout <- v_sf_out + v_sz_out + vep
    
x11()
layout(matrix(1:3,1,3))
matplot(cbind(q_sf_in,q_sf_out),type="l",main=paste("q_sf",nstep),xlim=c(0,200))
matplot(cbind(q_sz_in,q_sz_out),type="l",main=paste("q_sz",nstep),xlim=c(0,200))#)
plot(s[-n]+vin[-1]-vout[-1]-s[-1],xlim=c(0,200))#)
#matplot(cbind(s_sz,s_uz,s_rz,s_sf),type="l")

#plot(s[-n]+vin[-1]-vout[-1]-s[-1])
x11(); layout(matrix(1:4,2));
plot(s_sz,xlim=c(0,200));plot(s_uz,xlim=c(0,200));
plot(s_rz,xlim=c(0,200));plot(s_sf,xlim=c(0,200))


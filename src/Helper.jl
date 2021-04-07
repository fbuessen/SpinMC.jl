using MPI

function MPISendrecvFloat(sendbuf::Float64, dest::Integer, comm::MPI.Comm)::Float64
    sbuf = [sendbuf]
    rbuf = [sendbuf]
    MPI.Sendrecv!(sbuf, dest, 0, rbuf, dest, 0, comm)
    return rbuf[1]
end

function MPISendBool(sendbuf::Bool, dest::Integer, comm::MPI.Comm)
    sbuf = [sendbuf]
    MPI.Send(sbuf, dest, 0, comm)
    return nothing
end

function MPIRecvBool(dest::Integer, comm::MPI.Comm)::Bool
    rbuf = [ false ]
    MPI.Recv!(rbuf, dest, 0, comm)
    return rbuf[1]
end

function MPIBcastBool(sendbuf::Bool, root::Integer, comm::MPI.Comm)::Bool
    sbuf = [ sendbuf ]
    MPI.Bcast!(sbuf, root, comm)
    return sbuf[1]
end
{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 68,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 195,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "5\n",
      "[[[[0.  0.5 1.  1.5 2.  2.5 3.  3.5 4.  4.5]]]]\n"
     ]
    }
   ],
   "source": [
    "Lx, Ly, Lz = 5, 5, 5\n",
    "nodes = 5\n",
    "hx = Lx / nodes\n",
    "hy = Ly / nodes \n",
    "hz = Lz / nodes \n",
    "\n",
    "\n",
    "time_ticks = 10\n",
    "tau = 0.5\n",
    "\n",
    "X = np.arange(nodes)[:, np.newaxis, np.newaxis, np.newaxis] / (nodes - 1) * (Lx - 1)\n",
    "print(Lx)\n",
    "Y = np.arange(nodes)[np.newaxis, :, np.newaxis, np.newaxis] / (nodes - 1) * (Ly - 1)\n",
    "Z = np.arange(nodes)[np.newaxis, np.newaxis, :, np.newaxis] / (nodes - 1) * (Lz - 1)\n",
    "T = np.arange(time_ticks)[np.newaxis, np.newaxis, np.newaxis, :] * tau\n",
    "\n",
    "print(T)\n",
    "\n",
    "analytical = np.sin(3 * np.pi * X / Lx) * np.sin( 3 * np.pi * Y / Ly) * np.sin(np.pi * Z / Lz) * np.cos(4 * (T + np.pi)) \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 206,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1.0"
      ]
     },
     "execution_count": 206,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "hz"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 197,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[[[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [-0.    -0.    -0.    -0.    -0.   ]]\n",
      "\n",
      " [[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]]\n",
      "\n",
      " [[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]]]\n"
     ]
    }
   ],
   "source": [
    "np.set_printoptions(precision=3)\n",
    "\n",
    "print(analytical[:,:,:,1])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 198,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "0.0\n",
      "3.141592653589793 0 5 \n",
      "\n",
      "0.9510565162951536\n",
      "3.141592653589793 1 5 \n",
      "\n",
      "-0.587785252292473\n",
      "3.141592653589793 2 5 \n",
      "\n",
      "-0.5877852522924734\n",
      "3.141592653589793 3 5 \n",
      "\n",
      "0.9510565162951535\n",
      "3.141592653589793 4 5 \n",
      "\n"
     ]
    }
   ],
   "source": [
    "for i in range(5):\n",
    "    print(np.sin(3 * np.pi * i / Lx))\n",
    "    print(np.pi, i, Lx, '\\n')\n",
    "    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 211,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[[[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [-0.    -0.    -0.    -0.    -0.   ]]\n",
      "\n",
      " [[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [-0.    -0.085 -0.137 -0.137 -0.085]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]]\n",
      "\n",
      " [[-0.    -0.    -0.    -0.    -0.   ]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [ 0.     0.137  0.221  0.221  0.137]\n",
      "  [-0.    -0.221 -0.358 -0.358 -0.221]]]\n",
      "[[[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     1.243  2.011  2.011  0.   ]\n",
      "  [ 0.    -0.768 -1.243 -1.243  0.   ]\n",
      "  [ 0.    -0.768 -1.243 -1.243  0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.    -0.768 -1.243 -1.243  0.   ]\n",
      "  [ 0.     0.475  0.768  0.768  0.   ]\n",
      "  [ 0.     0.475  0.768  0.768  0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.    -0.768 -1.243 -1.243  0.   ]\n",
      "  [ 0.     0.475  0.768  0.768  0.   ]\n",
      "  [ 0.     0.475  0.768  0.768  0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]]\n",
      "\n",
      " [[ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]\n",
      "  [ 0.     0.     0.     0.     0.   ]]]\n"
     ]
    }
   ],
   "source": [
    "def laplasian(a):\n",
    "#     print(a)\n",
    "#     print('\\n\\n\\n\\n\\n\\n')\n",
    "    padded = np.pad(a, 1, 'constant')\n",
    "    out = (padded[:-2,1:-1,1:-1] - 2 * padded[1:-1,1:-1,1:-1] + padded[2:,1:-1,1:-1]) / hx ** 2 + \\\n",
    "          (padded[1:-1,:-2,1:-1] - 2 * padded[1:-1,1:-1,1:-1] + padded[1:-1,2:,1:-1]) / hy ** 2 + \\\n",
    "          (padded[1:-1,1:-1,:-2] - 2 * padded[1:-1,1:-1,1:-1] + padded[1:-1,1:-1,2:]) / hz ** 2\n",
    "    \n",
    "    out[0,:,:] = 0\n",
    "    out[-1,:,:] = 0\n",
    "    out[:,0,:] = 0\n",
    "    out[:,-1,:] = 0\n",
    "    out[:,:,0] = 0\n",
    "    out[:,:,-1] = 0\n",
    "    \n",
    "    return out\n",
    "\n",
    "result = np.zeros((nodes, nodes, nodes, time_ticks))\n",
    "result[:,:,:,0] = analytical[:,:,:,0]\n",
    "result[:,:,:,1] = analytical[:,:,:,1]\n",
    "\n",
    "\n",
    "# result[:,:,:,1] = result[:,:,:,0] + tau ** 2 * laplasian(result[:,:,:,0]) / 2\n",
    "# print(analytical[:,:,:,1])\n",
    "# print(analytical[:,:,:,1])\n",
    "\n",
    "\n",
    "for i in range(2, time_ticks):\n",
    "    print(result[:,:,:,i-1])\n",
    "    lpls = laplasian(result[:,:,:,i-1])\n",
    "#     result[:,:,:,i] = 2 * result[:,:,:,i-1] - result[:,:,:,i-2] + tau**2 * lpls\n",
    "    print(lpls)\n",
    "    break\n",
    "    print(np.max(np.abs(analytical[:,:,:,i] - result[:,:,:,i])))\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 164,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "92"
      ]
     },
     "execution_count": 164,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "np.max(np.array([[[0,1],[0,2]],[[1,2],[3,92]]]))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}

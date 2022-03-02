'''
科学计算8
韩秋原 
202111069028
'''
import math
import numpy as np
import numpy.ma as ma
import scipy.signal as signal
from matplotlib import pyplot as plt
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource

# from osgeo import gdal_array
class DEM:
    def __init__(self,data,info):
        self.raw_data = data
        self.size = info['cellsize']
        self.col = info['ncols']
        self.row = info['nrows']
        self.nan = info['NODATA_value']
        self.__xllcorner = info['xllcorner']
        self.__yllcorner = info['yllcorner']
        self.basic ={}
        self.window=[]
        #矩阵初始化
        self.__generate_basic_info()
    #获取统计数据    
    def __generate_basic_info(self):
        self.__masked_data = np.ma.masked_array(self.raw_data,mask=[self.raw_data==self.nan])
        self.__get_max(self.__masked_data)
        self.__get_min(self.__masked_data)
        self.__get_mean(self.__masked_data)
        self.__get_std(self.__masked_data)
        self.__figure_slope_numpy()
        self.__sorted_data=self.__sort_data(self.__masked_data)
        self.__sorted_slope=self.__sort_data(self.__masked_slope)
        self.__sorted_aspect=self.__sort_data(self.__masked_aspect)
        self.__write_asc()
        
    def __get_max(self,data):
        self.basic['max'] = np.max(data)
    def export_max(self):
        return self.basic['max']
    def __get_min(self,data):
        self.basic['min'] = np.min(data)
    def export_min(self):
        return self.basic['min']
    def __get_mean(self,data):
        self.basic['mean'] = np.average(data)
    def __get_std(self,data):
        self.basic['std'] = np.std(data)
    #绘制直方图
    def plot_hist(self,data_name='DEM',sub=[1,1,1]):
        if data_name=='DEM':
            y = self.__sorted_data
        elif data_name=='Slope':
            y = self.__sorted_slope
        elif data_name=='Aspect':
            y = self.__sorted_aspect
        a=4
        max_y=max(y)
        min_y=min(y)
        bins=10
#         n,bins,patches=plt.hist(y,bins)
        x=np.linspace(min_y,max_y,bins)
        return x,y,bins,max_y,min_y
    #高程分布画图
    def plot_data(self):
        #按数据标号绘制每个点
        y = self.__sorted_data
#         x = np.arange(y.shape[0])
        plt.xlabel("x axis")
        plt.ylabel("height")
        plt.plot(x,y)
        plt.show()
    #绘制矩阵
    def plot_data_2d(self,name='DEM'):
        if name=='DEM':
            data=self.__masked_data
        elif name=='Slope':
            data=self.__masked_slope
        elif name=='Aspect':
            data=self.__masked_aspect
#         plt.matshow(data)
#         plt.title('2D figure of '+name)
#         plt.xlabel("x axis")
#         plt.ylabel("y axis")
#         plt.colorbar()
#         plt.savefig("2D figure "+name)
        return data
    def plot_DEM_3D(self):
        pass
    #高程排序    
    def __sort_data(self,data):
#         data = data.reshape(1,-1)
#shape时（行数，列数）
#         data = data.flatten()
        data = data[~data.mask]
#         print(data.shape)
        return np.sort(data)

    def __aspect_degree_numpy(self):
        '''
        使用numpy计算角度
        '''
        self.__aspect_degrees=self.__slope_sn/self.__slope_we
        self.__aspect_degrees[(self.__slope_sn>0) & (self.__slope_we>0)]=np.degrees(np.arctan(self.__aspect_degrees[(self.__slope_sn>0) & (self.__slope_we>0)]))
        self.__aspect_degrees[(self.__slope_sn<=0) & (self.__slope_we>0)]=90+np.degrees(np.arctan(-1*self.__aspect_degrees[(self.__slope_sn<=0) & (self.__slope_we>0)]))
        self.__aspect_degrees[(self.__slope_sn>=0) & (self.__slope_we<0)]=270+np.degrees(np.arctan(-1*self.__aspect_degrees[(self.__slope_sn>=0) & (self.__slope_we<0)]))
        self.__aspect_degrees[(self.__slope_sn<0) & (self.__slope_we<0)]=180+np.degrees(np.arctan(self.__aspect_degrees[(self.__slope_sn<0) & (self.__slope_we<0)]))
        self.__aspect_degrees[(self.__slope_sn<0) & (self.__slope_we==0)]=180
        self.__aspect_degrees[(self.__slope_sn>0) & (self.__slope_we==0)]=0
#         return self.__aspect_degrees
    def __data_pre_process_2(self):
        data = self.raw_data
        #围2圈
        data = np.pad(data,((2,2),(2,2)),'constant',constant_values = self.nan)
        #拆矩阵
        for row in range(3):
            for col in range(3):
                arr = data[row:(row+data.shape[0]-2),col:(col+data.shape[1]-2)]
                self.window.append(np.ma.masked_array(arr,mask=[arr==self.nan]))     
        self.__e5 = self.window[0]
        self.__e2 = self.window[1]
        self.__e6 = self.window[2]
        self.__e1 = self.window[3]
        self.__e3 = self.window[5]
        self.__e8 = self.window[6]
        self.__e4 = self.window[7]
        self.__e7 = self.window[8]
    def __data_pre_process_1(self):
        '''
        数据预处理
        '''
        data = self.raw_data
        #矩阵偏移
        e3 = np.pad(data[:,2:],((1,1),(1,3)),'constant',constant_values = self.nan)
        e1 = np.pad(data[:,:],((1,1),(1,1)),'constant',constant_values = self.nan)
        e4 = np.pad(data[:,1:],((2,0),(1,2)),'constant',constant_values = self.nan)
        e2 = np.pad(data[1:,1:],((1,2),(1,2)),'constant',constant_values = self.nan)
        e5 = np.pad(data[:,:],((0,2),(1,1)),'constant',constant_values = self.nan)
        e6 = np.pad(data[:,2:],((0,2),(1,3)),'constant',constant_values = self.nan)
        e7 = np.pad(data[:,2:],((2,0),(1,3)),'constant',constant_values = self.nan)
        e8 = np.pad(data[:,:],((2,0),(1,1)),'constant',constant_values = self.nan)
        #掩膜
        self.__e1 = np.ma.masked_array(e1,mask=[e1==self.nan])
        self.__e2 = np.ma.masked_array(e2,mask=[e2==self.nan])
        self.__e3 = np.ma.masked_array(e3,mask=[e3==self.nan])
        self.__e4 = np.ma.masked_array(e4,mask=[e4==self.nan])
        self.__e5 = np.ma.masked_array(e5,mask=[e5==self.nan])
        self.__e6 = np.ma.masked_array(e6,mask=[e6==self.nan])
        self.__e7 = np.ma.masked_array(e7,mask=[e7==self.nan])
        self.__e8 = np.ma.masked_array(e8,mask=[e8==self.nan])
    def figure_slope_scipy(self,arg=2):
        '''
        使用scipy卷积，速度更快，只做缺省值补全，不处理悬崖部分
        '''
        data = self.raw_data
        data = np.ma.masked_array(data,mask=[data==self.nan])
        kernel_1_we=[[0,0,0],[1,0,-1],[0,0,0]]
        kernel_1_sn=[[0,-1,0],[0,0,0],[0,1,0]]
        kernel_2_we=[[1,0,-1],[2,0,-2],[1,0,-1]]
        kernel_2_sn=[[-1,-2,-1],[0,0,0],[1,2,1]]
        if arg==1:
            self.__slope_we=signal.convolve2d(data,kernel_1_we)/(2*self.size)
            self.__slope_sn=signal.convolve2d(data,kernel_1_sn)/(2*self.size)
        elif arg==2:
            self.__slope_we=signal.convolve2d(data,kernel_2_we)/(8*self.size)
            self.__slope_sn=signal.convolve2d(data,kernel_2_sn)/(8*self.size)
        self.__calculation()
        self.__slope[(self.__slope==0) | (self.__slope>80)]=self.nan
        self.__aspect[(self.__aspect==0)|(self.__aspect>3)]=self.nan
            
    def __figure_slope_numpy(self,pre=1,arg=2,angle=True):
        '''
        使用numpy计算坡度坡向原始值
        '''
        if pre==1:
            self.__data_pre_process_1()
        elif pre==2:
            self.__data_pre_process_2()
        if arg==1:
            self.__slope_we = ((self.__e1-self.__e3)/(2*self.size))[1:-1,1:-1]
            self.__slope_sn = ((self.__e4-self.__e2)/(2*self.size))[1:-1,1:-1]
        elif arg==2:
            self.__slope_we = (((self.__e8+2*self.__e1+self.__e5)-(self.__e7+2*self.__e3+self.__e6))/(8*self.size))[1:-1,1:-1]
            self.__slope_sn = (((self.__e8+2*self.__e4+self.__e7)-(self.__e5+2*self.__e2+self.__e6))/(8*self.size))[1:-1,1:-1]
        self.__calculation(angle)

    def __calculation(self,angle=True):
        '''
        通过sn和we计算坡度和坡向
        同时处理scipy的ndarray
        numpy的MaskedArray
        '''
        self.__slope = np.degrees(np.arctan(np.sqrt(np.power(self.__slope_sn,2)+np.power(self.__slope_we,2))))
        if angle:
            self.__aspect_degree_numpy()
            #暂定解决方案1（因为写入没给degrees选项）
            self.__aspect=self.__aspect_degrees
        else:
            self.__aspect = np.arctan2(self.__slope_we,self.__slope_sn)
            #mask前的计算，we可能出现0
            #废弃的方法
#         self.__aspect = np.divide(self.__slope_sn,self.__slope_we,out = np.zeros_like(self.__slope_sn),where=self.__slope_we!=0)
        if type(self.__slope) is np.ma.core.MaskedArray:
        #保存mask数据和恢复后的数据
            self.__masked_slope = self.__slope
            self.__masked_aspect = self.__aspect
            self.__slope = self.__slope.filled(fill_value=self.nan)
            self.__aspect = self.__aspect.filled(fill_value=self.nan)
    #写文件
    def __write_asc(self):
        '''
        文件名slope,aspect
        '''
        ncols = self.__slope.shape[0]
        nrows = self.__slope.shape[1]
        head = 'ncols\t'+str(ncols)+'\n'+\
        'nrows\t'+str(nrows)+'\n'+\
        'xllcorner\t'+str(self.__xllcorner)+'\n'+\
        'yllcorner\t'+str(self.__yllcorner)+'\n'+\
        'cellsize\t'+str(self.size)+'\n'+\
        'NODATA_value\t'+str(self.nan)+'\n'
        with open('slope.asc','w') as f:
            f.write(head)
            np.savetxt(f,self.__slope,delimiter='\t',fmt = '%.02f')
        with open('aspect.asc','w') as f:
            f.write(head)
            np.savetxt(f,self.__aspect,delimiter='\t',fmt = '%.02f')
# 读文件
def read_asc(name):
    with open(name,'r') as f:
        info =dict()
        info['ncols']=int(f.readline().split()[1])
        info['nrows']=int(f.readline().split()[1])
        info['xllcorner']=float(f.readline().split()[1])
        info['yllcorner']=float(f.readline().split()[1])
        info['cellsize']=int(f.readline().split()[1])
        info['NODATA_value']=int(f.readline().split()[1])
    #读入矩阵
    data = np.loadtxt(name,skiprows=6)
    return info,data
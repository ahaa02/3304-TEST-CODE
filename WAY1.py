from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
import math
import time

class Window2d():
    model2d = None

    def __init__(self, m):
        self.model2d = m

    def update(self):
        self.draw2d()

    def drawModel2d(self): # 출력 형식 초기 설정
       
        for f in self.model2d.faces:
            glColor3f(f.colorR, f.colorG, f.colorB)
            he = f.halfedge
            count = 0

            glBegin(GL_LINE_LOOP)
            while True:
                glVertex3f(he.vertex.x, he.vertex.y, he.vertex.z)
                he = he.next
                count += 1
                if count == f.NumOfEdges:
                    break
            glEnd()


    def draw2d(self): # 출력 형식 설정
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        ##set camera
        maxX = -1000000
        minX = 1000000
        maxY = -1000000
        minY = 1000000
        maxZ = -1000000
        minZ = 1000000
        for v in self.model2d.vertices:
            if v.x > maxX:
                maxX = v.x
            if v.x < minX:
                minX = v.x
            if v.y > maxY:
                maxY = v.y
            if v.y < minY:
                minY = v.y
            if v.z > maxZ:
                maxZ = v.z
            if v.z < minZ:
                minZ = v.z
        gluLookAt((minX+maxX)/2, (minY+maxY)/2, 1.3*max((maxY-minY), (maxX-minX)) + maxZ, (minX+maxX)/2, (minY+maxY)/2, (minZ+maxZ)/2, 0.0, 1.0, 0.0)
        self.drawModel2d()
        glFlush() 

    def init2d(self, width, height):
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glEnable(GL_DEPTH_TEST)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(width)/float(height), 0.1, 10000.0)

    def resize2d(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(w)/float(h), 0.1, 10000.0)

    def mouseInput(self, button, state, x, y):
        self.update()

    def keyInput(self, key, x, y):
        if key == b' ':  # 스페이스 입력 시 폴리곤 상태 출력
            self.polygonsInfo()
        if key == b'p': # p 입력 시 평면에 폴리곤 출력
            self.toPlane()
        if key == b's': # s 입력 시 세타 값을 고려한 폴리곤 출력
            self.calucSita()
        self.update()

    def calucW(self, sitas, idVtm, idEtm): 
        np.set_printoptions(suppress=True, threshold = 1000000, linewidth=200)
        c = np.zeros((idVtm*2, idEtm))
        ws = np.zeros((idEtm, 1))
        b = np.zeros((idVtm*2, 1))

        for vtm in self.model2d.vtms:
            sita_sum = 0
            for i in range(len(vtm.halfedges)):
                he = vtm.halfedges[i]
                if (i % 2 == 1):
                    sita_sum += he.vertex_next.sita_large
                    if he.pair == None:
                        c[vtm.id*2, he.id] = math.cos(sita_sum)
                        c[vtm.id*2+1, he.id] = math.sin(sita_sum)
                    else:
                        if vtm.id < he.vertex.vtm.id:
                            c[vtm.id*2, he.id] = math.cos(sita_sum)
                            c[vtm.id*2+1, he.id] = math.sin(sita_sum)
                        else:
                            c[vtm.id*2, he.id] = math.cos(sita_sum)
                            c[vtm.id*2+1, he.id] = math.sin(sita_sum)
                            b[vtm.id*2, 0] += math.cos(sita_sum)*2*self.getLength(he.vertex, he.vertex_next)*math.sin(-0.5*sitas[he.id, 0])
                            b[vtm.id*2+1, 0] += math.sin(sita_sum)*2*self.getLength(he.vertex, he.vertex_next)*math.sin(-0.5*sitas[he.id, 0])

        cPlus = np.linalg.pinv(c)
        ws = np.dot(cPlus, b) + np.dot((np.eye(c.shape[1]) - np.dot(cPlus, c)), ws)
        print("norm(C*ws - b) = " + str(np.linalg.norm(np.dot(c,ws)-b)))


    def calucSita(self):  # 세타 값 계산 함수

        np.set_printoptions(suppress=True, threshold = 1000000, linewidth=200)
        idEtm = 0
        for he in self.model2d.halfedges:
            if he.id == None:
                if he.pair == None:
                    he.id = idEtm
                else:
                    he.id = idEtm
                    he.pair.id = idEtm
                idEtm += 1
        idVtm = 0
        for vtm in self.model2d.vtms:
            vtm.id = idVtm
            idVtm += 1

        c = np.zeros((idVtm, idEtm))
        sitas = np.zeros((idEtm, 1))
        g = np.zeros((idVtm, 1))

        for vtm in self.model2d.vtms:
            gValue = 2*math.pi
            for i in range(len(vtm.halfedges)):
                he = vtm.halfedges[i]
                if (i % 2 == 1):
                    if he.pair == None:
                        c[vtm.id, he.id] = 1
                    else:
                        if (vtm.id < he.vertex.vtm.id):
                            c[vtm.id, he.id] = 1
                        else:
                            c[vtm.id, he.id] = -1
                    sita = self.getSita(he.next.vertex, he.next.vertex_next, he.vertex_next, he.vertex)
                    he.vertex_next.sita = sita
                    gValue -= sita
            g[vtm.id, 0] = gValue

        cPlus = np.linalg.pinv(c)
        IminusCC = np.eye(c.shape[1]) - np.dot(cPlus, c)
        sitas = np.dot(cPlus, g) + np.dot(IminusCC, sitas)

        iter = 0
        while True:
            E = 0
            deltaE = np.zeros((idEtm, 1))
            for vtm in self.model2d.vtms:
                sita_large_right = None
                idx_sita_large_right = None
                flag_rightEtmRev = False
                if vtm.halfedges[0].pair == None:
                    sita_large_right = sitas[vtm.halfedges[-1].id, 0]
                    flag_rightEtmRev = False
                    idx_sita_large_right = vtm.halfedges[-1].id
                else:
                    if (vtm.id < vtm.halfedges[0].vertex_next.vtm.id):
                        sita_large_right = sitas[vtm.halfedges[0].id, 0]
                        flag_rightEtmRev = False
                        idx_sita_large_right = vtm.halfedges[0].id
                    else:
                        sita_large_right = -1 * sitas[vtm.halfedges[0].id, 0]
                        flag_rightEtmRev = True
                        idx_sita_large_right = vtm.halfedges[0].id

                for i in range(len(vtm.halfedges)):
                    he = vtm.halfedges[i]
                    sita_etm = None
                    sita_boarder = None
                    sita_large = None
                    sita_large_left = None
                    flag_leftEtmRev = False
                    if (i % 2 == 1):
                        if he.pair == None:
                            sita_boarder = sitas[he.id, 0]
                            sita_large_left = sita_boarder
                            flag_leftEtmRev = False
                        else:
                            if (vtm.id < he.vertex.vtm.id):
                                sita_etm = sitas[he.id, 0]
                                sita_large_left = sita_etm
                                flag_leftEtmRev = False
                            else:
                                sita_etm = -1 * sitas[he.id, 0]
                                sita_large_left = sita_etm
                                flag_leftEtmRev = True

                        sita_large = sita_large_right/2 + he.vertex_next.sita + sita_large_left/2
                        he.vertex_next.sita_large = sita_large

                        k = math.pi

                        if sita_etm != None:
                            if sita_etm <= -1 * math.pi:
 
                                E += (sita_etm + math.pi)**2
                                if flag_leftEtmRev:
                                    deltaE[he.id, 0] = -1 * (2*sita_etm + 2*math.pi - 2*k)
                                else:
                                    deltaE[he.id, 0] = 2*sita_etm + 2*math.pi - 2*k
                            elif sita_etm >= math.pi:
 
                                E += (sita_etm - math.pi)**2
                                if flag_leftEtmRev:
                                    deltaE[he.id, 0] = -1 * (2*sita_etm - 2*math.pi + 2*k)
                                else:
                                    deltaE[he.id, 0] = 2*sita_etm - 2*math.pi + 2*k
                            else:
                                x = 1
                        if sita_large != None and he.pair != None and vtm.halfedges[0].pair != None:
                            if sita_large < 0:
  
                                E += (sita_large)**2
                                if flag_leftEtmRev:
                                    deltaE[he.id, 0] = -1 * (0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - 0.5*k)
                                else:
                                    deltaE[he.id, 0] = 0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - 0.5*k
                                if flag_rightEtmRev:
                                    deltaE[idx_sita_large_right, 0] = -1 * (0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - 0.5*k)
                                else:
                                    deltaE[idx_sita_large_right, 0] = 0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - 0.5*k
                            elif sita_large >= math.pi:
                                E += (sita_large - math.pi)**2
                                if flag_leftEtmRev:
                                    deltaE[he.id, 0] = -1 * (0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - math.pi + 0.5*k)
                                else:
                                    deltaE[he.id, 0] = 0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - math.pi + 0.5*k
                                if flag_rightEtmRev:
                                    deltaE[idx_sita_large_right, 0] = -1 * (0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - math.pi + 0.5*k)
                                else:
                                    deltaE[idx_sita_large_right, 0] = 0.5*sita_large_right + 0.5*sita_large_left + he.vertex_next.sita - math.pi + 0.5*k
                            else:
                                x = 1
                        if sita_boarder != None:
                            if sita_boarder < math.pi:
                                E += (sita_boarder - math.pi)**2
                                if flag_leftEtmRev:
                                    deltaE[he.id, 0] = -1 * (2*sita_boarder - 2*math.pi - 2*k)
                                else:
                                    deltaE[he.id, 0] = 2*sita_boarder - 2*math.pi - 2*k
                            else:
                                x = 1
                        idx_sita_large_right = he.id
                        sita_large_right = sita_large_left
                        flag_rightEtmRev = flag_leftEtmRev
            iter += 1
            alpha = 0.1
            print("E_sitas = " + str(E))
            if E < 0.00000000000001 or iter > 100:
                break
            else:
                deltaE = np.dot(IminusCC, deltaE)
                sitas = sitas - alpha*deltaE
                self.update()

        self.setBySita(self.model2d.faces[5], sitas)
        self.update()
        self.calucW(sitas, idVtm, idEtm)


    def setBySita(self, f, sitas):
        f.flag_setBySita = True
        heStart = f.halfedge
        he = heStart
        while True:
            if he.pair != None:
                if not he.pair.face.flag_setBySita:
                    if he.vertex.vtm.id < he.vertex_next.vtm.id:
                        sita_will = -1 * sitas[he.id, 0]
                    elif he.vertex.vtm.id > he.vertex_next.vtm.id:
                        sita_will = sitas[he.id, 0]
                    else:
                   
                    sita_now = self.getSita(he.vertex_next, he.vertex, he.pair.vertex, he.pair.vertex_next)
                    sita_rotate = sita_will - sita_now
                  
                    a = he.pair.vertex_next.x
                    b = he.pair.vertex_next.y
                    c = he.pair.vertex.x
                    d = he.pair.vertex.y
                    x = c + (a-c)*math.cos(sita_rotate) - (b-d)*math.sin(sita_rotate)
                    y = d + (a-c)*math.sin(sita_rotate) + (b-d)*math.cos(sita_rotate)
                    he.pair.vertex_next.x = x
                    he.pair.vertex_next.y = y
               
                    a = he.pair.next.vertex_next.x
                    b = he.pair.next.vertex_next.y
                    c = he.pair.vertex.x
                    d = he.pair.vertex.y
                    x = c + (a-c)*math.cos(sita_rotate) - (b-d)*math.sin(sita_rotate)
                    y = d + (a-c)*math.sin(sita_rotate) + (b-d)*math.cos(sita_rotate)
                    he.pair.next.vertex_next.x = x
                    he.pair.next.vertex_next.y = y
            
                    self.setBySita(he.pair.face, sitas)
            he = he.next
            if he == heStart:
                break


    def polygonsInfo(self):
        print("----------model2d-----------")
        print("len(faces) = " + str(len(self.model2d.faces)))
        print("len(halfedges) = " + str(len(self.model2d.halfedges)))
        print("len(vertices) = " + str(len(self.model2d.vertices)))
        print("len(vtms) = " + str(len(self.model2d.vtms)))

    def toPlane(self):
        for f in self.model2d.faces:
            if f.NumOfEdges == 3:
                he = f.halfedge
                length_next = self.getLength(he.vertex, he.vertex_next)
                length_prev = self.getLength(he.prev.vertex_next, he.prev.vertex)
                sita = self.getSita(he.vertex, he.vertex_next, he.prev.vertex_next, he.prev.vertex)
                if sita > math.pi:
                    sita = 2*math.pi - sita
                he.vertex.x = 0
                he.vertex.y = 0
                he.vertex.z = 0
                he.vertex_next.x = length_next
                he.vertex_next.y = 0
                he.vertex_next.z = 0
                he.prev.vertex.x = length_prev * math.cos(sita)
                he.prev.vertex.y = length_prev * math.sin(sita)
                he.prev.vertex.z = 0
            else:
                print("face.NumOfEdges != 3")


    def getSita(self, v1, v2, v3, v4):
        v1v2_length = self.getLength(v1,v2)
        v3v4_length = self.getLength(v3,v4)
        sitaCos = ((v2.x - v1.x) * (v4.x - v3.x) + (v2.y - v1.y) * (v4.y - v3.y) + + (v2.z - v1.z) * (v4.z - v3.z)) / (v1v2_length * v3v4_length)
        sita = math.acos(sitaCos)
        if self.getGaisekiZ(v1, v2, v3, v4) < 0:
            sita = 2*math.pi - sita
        return sita

    def getLength(self, v1, v2):
        return math.sqrt(((v2.x - v1.x) ** 2) + ((v2.y - v1.y) ** 2) + ((v2.z - v1.z) ** 2))

    def getGaisekiZ(self, v1, v2, v3, v4):

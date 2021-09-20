const { create } = require('domain');
const {app, BrowserWindow, ipcMain} = require('electron');
const fs = require('fs');
const path = require('path');
const {spawn} = require('child_process');
const isDev=require('electron-is-dev');


let mainWindow;

const createWindow=()=>
{
    mainWindow = new BrowserWindow(
            {
                width:1440,
                height:900,
                webPreferences:
                {
                    nodeIntegration:false,
                    contextIsolation:true,
                    enableRemoteModule:false,
                    nativeWindowOpen:true,
                    preload:path.join(__dirname,"preload.js")
                },
            }
        );
    // mainWindow.webContents.openDevTools();
    mainWindow.loadFile(path.join(__dirname,"preload.js"));
    const startUrl= isDev ? 'http://localhost:3000' : `file://${path.join(__dirname,'../build/index.html')}`
    mainWindow.loadURL('http://localhost:3000');

}


app.on('ready',createWindow);


const readPlot = (arg)=>
{
    const py = spawn('python',[path.join(__dirname,"/python_code/3_plot_star_interactions.py"),arg]);
    py.stdout.on('close',()=>
    {
        fs.readFile(path.join(__dirname,'/python_code/VelocityOrionStars2.jpg'),(err,data)=>
        {
            if(err)
                console.log(err);
            else
                {
                    mainWindow.webContents.send('send_plot',data);
                }
        });
    });
}

ipcMain.on('get_plot',(event,arg)=>
{   
    readPlot(arg);
})


